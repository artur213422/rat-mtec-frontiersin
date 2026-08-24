from __future__ import annotations

from dataclasses import dataclass

import polars as pl


@dataclass(frozen=True, slots=True)
class PreparedOrderedData:
    categories: list[str]
    y_values: list[float]


def prepare_data(picked_genes: list[str], tec_df: pl.DataFrame, with_orthologs: bool) -> PreparedOrderedData:
    required_columns: set[str] = {
        "gene_symbol",
        "log2_fold_change",
        "ctrl_expression_percentile",
    }

    if with_orthologs:
        required_columns.add("has_clean_ortholog")

    missing_columns: set[str] = required_columns - set(tec_df.columns)

    if missing_columns:
        raise ValueError(
            f"Missing required columns: {sorted(missing_columns)}"
        )

    if not picked_genes:
        return PreparedOrderedData(
            categories=[],
            y_values=[],
        )

    filter_expression: pl.Expr = (
        pl.col("gene_symbol").is_in(picked_genes)
    )

    if with_orthologs:
        filter_expression &= (
            pl.col("has_clean_ortholog").fill_null(False)
        )

    prepared_df: pl.DataFrame = (
        tec_df
        .filter(filter_expression)
        .filter(
            pl.col("gene_symbol").is_not_null()
            & pl.col("log2_fold_change").is_not_null()
            & pl.col("log2_fold_change").is_finite()
            & pl.col("ctrl_expression_percentile").is_not_null()
            & pl.col("ctrl_expression_percentile").is_finite()
        )
        # Strongest negative KO-vs-CTRL fold change first
        .sort(
            "log2_fold_change",
            descending=False,
        )
    )

    categories: list[str] = (
        prepared_df
        .get_column("gene_symbol")
        .cast(pl.String)
        .to_list()
    )

    y_values: list[float] = (
        prepared_df
        .get_column("ctrl_expression_percentile")
        .cast(pl.Float64)
        .to_list()
    )

    return PreparedOrderedData(
        categories=categories,
        y_values=y_values,
    )