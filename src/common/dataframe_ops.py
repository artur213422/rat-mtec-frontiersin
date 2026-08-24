import re

import numpy as np
import polars as pl

from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.types.dataset import FloatMatrix
from src.common.types.enums import Genotypes
from src.common.utils.col_names import GENE_SYMBOL_COLUMN, ORTHOLOG_ID_COLUMN


def get_custom_genotype_expression_values(row: dict[str, object], pattern: str) -> list[float]:
    compiled_pattern: re.Pattern[str] = re.compile(
        rf"^{re.escape(pattern)}\d+$",
        flags=re.IGNORECASE,
    )

    return [
        float(value)
        for col, value in row.items()
        if compiled_pattern.fullmatch(col) and value is not None
    ]

def get_genotype_expression_values(row: dict[str, object], gt: Genotypes,) -> list[float]:
    pattern = re.compile(rf"^{gt.value}\d+$")

    return [
        float(value)
        for col, value in row.items()
        if pattern.fullmatch(col) and value is not None
    ]


def get_genotype_col_names(tec_df: pl.DataFrame, gt: Genotypes,) -> list[str]:
    prefix: str = gt.value.lower()

    matching_columns: list[str] = [
        column
        for column in tec_df.columns
        if column.lower().startswith(prefix)
    ]

    if not matching_columns:
        raise ValueError(
            f"No columns found for genotype {gt!r}. "
            f"Expected columns starting with {prefix!r}. "
            f"Available columns: {tec_df.columns}"
        )

    return matching_columns






def print_null_gene_summary(data: pl.DataFrame, id_columns: tuple[str, ...]) -> None:
    sample_cols: list[str] = [
        column
        for column in data.columns
        if column not in id_columns
        and column != GENE_SYMBOL_COLUMN
        and column != ORTHOLOG_ID_COLUMN
        and not column.endswith("_gene_symbol")
        and not column.endswith("_ortholog_id")
    ]

    if not sample_cols:
        print(f"Total genes: {data.height}")
        print("No sample columns found")
        return

    non_null_value_count: pl.Expr = pl.sum_horizontal(
        pl.col(sample_cols).is_not_null().cast(pl.Int64)
    )

    genes_with_nulls: int = data.filter(
        pl.any_horizontal(pl.col(sample_cols).is_null())
    ).height

    fully_populated_genes: int = data.filter(
        pl.all_horizontal(pl.col(sample_cols).is_not_null())
    ).height

    genes_with_less_than_5_values: int = data.filter(
        non_null_value_count < 5
    ).height

    total_genes: int = data.height

    print(f"Total genes: {total_genes}")
    print(f"Genes containing nulls: {genes_with_nulls}")
    print(f"Fully populated genes: {fully_populated_genes}")
    print(f"Genes with less than 5 values: {genes_with_less_than_5_values}")


def to_log2(data: pl.DataFrame, selected_col: str) -> pl.DataFrame:
    return data.with_columns(
        pl.col(selected_col)
        .log(base=2)
        .alias(selected_col)
    )


def to_linear(data: pl.DataFrame, selected_col: str) -> pl.DataFrame:
    return data.with_columns(
        pl.lit(2.0)
        .pow(pl.col(selected_col))
        .alias(selected_col)
    )

def extract_matrix(data: pl.DataFrame, gt: Genotypes) -> FloatMatrix:
    rows: list[list[float]] = []

    for row in data.iter_rows(named=True):
        values: list[float] = []

        match gt:
            case Genotypes.KO:
                values = get_genotype_expression_values(row=row, gt=Genotypes.KO)

            case Genotypes.CTRL:
                values = get_genotype_expression_values(row=row, gt=Genotypes.CTRL)

            case Genotypes.ALL:
                ctrl_values: list[float] = get_genotype_expression_values(row=row, gt=Genotypes.CTRL)

                ko_values = get_genotype_expression_values(row=row, gt=Genotypes.KO)

                values = [
                    *ctrl_values,
                    *ko_values,
                ]

            case _:
                raise ValueError(f"Unsupported genotype: {gt}")

        rows.append(values)

    if not rows:
        raise ValueError("Cannot extract a matrix from an empty DataFrame.")

    expected_sample_count: int = len(rows[0])

    if expected_sample_count == 0:
        raise ValueError(f"No expression columns found for genotype {gt}.")

    for row_index, values in enumerate(rows):
        if len(values) != expected_sample_count:
            raise ValueError("All genes must have the same number of expression values. Row 0 has {expected_sample_count}, but row {row_index} has {len(values)}.")

    matrix: FloatMatrix = np.asarray(rows, dtype=np.float64)


    return matrix