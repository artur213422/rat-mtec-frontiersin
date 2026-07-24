import csv
from io import StringIO
from pathlib import Path

import polars as pl

from src.common.file_io import _open_text


def read_gsm_file(file: Path) -> pl.DataFrame:
    table_lines: list[str] = []
    in_table = False

    with _open_text(file) as handle:
        for line in handle:
            stripped = line.strip()

            if stripped == "!series_matrix_table_begin":
                in_table = True
                continue

            if stripped == "!series_matrix_table_end":
                break

            if in_table:
                table_lines.append(line)

    if not table_lines:
        raise ValueError(
            f"No GEO series matrix table found in file: {file}. "
            "Expected block between !series_matrix_table_begin and "
            "!series_matrix_table_end."
        )

    df = pl.read_csv(
        StringIO("".join(table_lines)),
        separator="\t",
        has_header=True,
        quote_char='"',
        null_values=["", "NA", "NaN", "null"],
        infer_schema_length=1_000,
        schema_overrides={"ID_REF": pl.Utf8},
    )

    expression_columns = [col for col in df.columns if col != "ID_REF"]

    return df.with_columns(
        pl.col(expression_columns).cast(pl.Float64, strict=False)
    )


def get_annotation_table(annot_file: Path) -> pl.DataFrame:
    annot = pl.read_csv(
        annot_file,
        separator="\t",
        comment_prefix="#",
        infer_schema_length=0,
    )

    if annot.is_empty():
        raise ValueError(f"Annotation table is empty: {annot_file}")

    return annot
