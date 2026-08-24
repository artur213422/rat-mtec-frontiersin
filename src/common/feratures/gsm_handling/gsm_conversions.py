import csv
import re
from collections import defaultdict
from dataclasses import dataclass

import polars as pl

from pathlib import Path

from src.common.feratures.gsm_handling.gsm_io import get_annotation_table, read_gsm_file
from src.common.paths import find_file_in_project


def get_sample_names(file_path: Path) -> list[str]:
    with file_path.open(encoding="utf-8") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                row = next(csv.reader([line], delimiter="\t", quotechar='"'))
                return row[1:]

    raise ValueError("!Sample_title not found.")


def get_platform_id(file_path: Path) -> str:
    with file_path.open(encoding="utf-8") as f:
        for line in f:
            if line.startswith("!Series_platform_id"):
                return line.rstrip().split("\t")[1].strip('"')

    raise ValueError("!Series_platform_id not found.")

def _get_gene_symbol(probe_id: str, annot: pl.DataFrame, selected_cols: list[str]) -> str | None:
    selected_col: str | None = next(
        (col for col in selected_cols if col in annot.columns),
        None,
    )

    if selected_col is None:
        raise ValueError(
            f"None of the requested columns found. "
            f"Requested: {selected_cols}. "
            f"Available: {annot.columns}"
        )

    result = (
        annot
        .filter(pl.col("ID") == probe_id)
        .select(selected_col)
    )


    if result.height == 0:
        return None

    gene_symbol: str | None = result.item()

    if gene_symbol is None or gene_symbol == "":
        return None

    return gene_symbol

@dataclass
class AnnotationResult:
    table: pl.DataFrame
    missing_gene_count: int = 0
    mapped_gene_count: int = 0

def annotate_expression_matrix(annotation_table_file: Path, expression_matrix: Path) -> AnnotationResult:

    try:
        matrix: pl.DataFrame = read_gsm_file(expression_matrix)
    except ValueError:
        matrix: pl.DataFrame = pl.read_csv(expression_matrix, separator="\t")

    print(f"Processing: {expression_matrix.stem}")

    annot_table: pl.DataFrame = get_annotation_table(annotation_table_file)

    gene_symbols: list[str | None] = []

    missing_gene_count: int = 0
    mapped_gene_count: int = 0

    for probe_id in matrix["ID_REF"]:

        selected_cols: list[str] = ["gene_assignment", "Gene Symbol"]
        gene_name: str | None = _get_gene_symbol(str(probe_id), annot_table, selected_cols)  ## Or Gene Symbol

        if gene_name is None:
            missing_gene_count += 1
            gene_symbols.append(None)
            continue

        mapped_gene_count += 1
        gene_symbols.append(gene_name)

    matrix = matrix.with_columns(
        pl.Series("gene_symbol", gene_symbols)
    )

    return AnnotationResult(
        table=matrix,
        missing_gene_count=missing_gene_count,
        mapped_gene_count=mapped_gene_count,
    )

def clean_sample_names(sample_names: list[str]) -> list[str]:
    """Remove the GEO suffix and append replicate numbers."""

    counts: defaultdict[str, int] = defaultdict(int)
    cleaned: list[str] = []

    for name in sample_names:
        base = re.sub(r"_4MJW\d+$", "", name)
        counts[base] += 1
        cleaned.append(f"{base}_{counts[base]}")

    return cleaned

def rename_expression_matrix_col(
    matrix: pl.DataFrame,
    col_names: list[str],
) -> pl.DataFrame:
    protected_cols: set[str] = {"ID_REF", "gene_symbol"}

    sample_cols: list[str] = [
        col for col in matrix.columns
        if col not in protected_cols
    ]

    if len(sample_cols) != len(col_names):
        raise ValueError(
            f"Number of sample columns and new names does not match. "
            f"Sample columns: {len(sample_cols)}, new names: {len(col_names)}"
        )

    rename_map: dict[str, str] = dict(zip(sample_cols, col_names, strict=True))

    matrix = matrix.rename(rename_map)

    # Put columns in a clean order.
    ordered_cols: list[str] = [
        "ID_REF",
        "gene_symbol",
        *col_names,
    ]

    missing_cols: list[str] = [
        col for col in ordered_cols
        if col not in matrix.columns
    ]

    if missing_cols:
        raise ValueError(f"Missing expected columns after rename: {missing_cols}")

    return matrix.select(ordered_cols)