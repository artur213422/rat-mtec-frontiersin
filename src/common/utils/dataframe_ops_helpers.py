import polars as pl
import re

from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.types.enums import Genotypes
from src.common.utils.col_names import GENE_SYMBOL_COLUMN, ORTHOLOG_ID_COLUMN


def _get_genotype_columns(df: pl.DataFrame, gt: Genotypes) -> list[str]:
    pattern: re.Pattern[str] = re.compile(rf"^{re.escape(gt.value)}\d+$")

    return [
        column
        for column in df.columns
        if pattern.fullmatch(column) is not None
    ]


def _get_clean_ortholog_filter(df: pl.DataFrame) -> pl.Expr:
    if "has_clean_ortholog" in df.columns:
        column: str = "has_clean_ortholog"
    elif "has_ortholog" in df.columns:
        column = "has_ortholog"
    else:
        raise ValueError("Missing clean ortholog column. Expected 'has_clean_ortholog' or 'has_ortholog'.")

    if df.schema[column] == pl.Boolean:
        return pl.col(column)

    return (
        pl.col(column)
        .cast(pl.Utf8)
        .str.to_lowercase()
        .is_in(["true", "1", "yes"])
    )

