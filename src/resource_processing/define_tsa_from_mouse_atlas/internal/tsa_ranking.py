from __future__ import annotations

import polars as pl


def rank_tsa_group(matrix: pl.DataFrame, group_name: str, tissue_groups: dict[str, list[str]],) -> pl.DataFrame:
    if "gene_symbol" not in matrix.columns:
        raise ValueError("Missing required column: gene_symbol")

    if group_name not in tissue_groups:
        raise ValueError("Unknown tissue group '{group_name}'. Available groups: {list(tissue_groups)}")

    group_tissues: list[str] = tissue_groups[group_name]

    if not group_tissues:
        raise ValueError(f"Tissue group '{group_name}' contains no tissues")

    missing_tissues: list[str] = [
        tissue
        for tissue in group_tissues
        if tissue not in matrix.columns
    ]

    if missing_tissues:
        raise ValueError(f"Tissue group '{group_name}' contains tissues missing from matrix: {missing_tissues}")

    tissue_cols: list[str] = [
        col
        for col in matrix.columns
        if col != "gene_symbol"
    ]

    group_tissue_set: set[str] = set(group_tissues)

    outside_tissues: list[str] = [
        tissue
        for tissue in tissue_cols
        if tissue not in group_tissue_set
    ]

    if not outside_tissues:
        raise ValueError(
            f"Tissue group '{group_name}' contains all available tissues"
        )

    group_exprs: list[pl.Expr] = [
        pl.col(tissue)
        for tissue in group_tissues
    ]

    outside_exprs: list[pl.Expr] = [
        pl.col(tissue)
        for tissue in outside_tissues
    ]

    ranked: pl.DataFrame = (
        matrix
        .with_columns(
            pl.max_horizontal(group_exprs).alias("group_max_log2"),
            pl.mean_horizontal(group_exprs).alias("group_mean_log2"),
            pl.max_horizontal(outside_exprs).alias("outside_max_log2"),
            pl.mean_horizontal(outside_exprs).alias("outside_mean_log2"),
        )
        .with_columns(
            (
                pl.col("group_max_log2")
                - pl.col("outside_max_log2")
            ).alias("max_specificity_log2"),

            (
                pl.col("group_mean_log2")
                - pl.col("outside_mean_log2")
            ).alias("mean_enrichment_log2"),

            pl.mean_horizontal(
                [
                    (
                        pl.col(tissue)
                        - pl.col("outside_max_log2")
                    ).clip(lower_bound=0.0)
                    for tissue in group_tissues
                ]
            ).alias("tsa_score"),
        )
        .sort("tsa_score", descending=True)
        .with_row_index("rank", offset=1)
        .select(
            "gene_symbol",
            pl.lit(group_name).alias("tissue_group"),
            "rank",
            "tsa_score",
            "group_max_log2",
            "group_mean_log2",
            "outside_max_log2",
            "max_specificity_log2",
            "mean_enrichment_log2",
        )
    )

    return ranked