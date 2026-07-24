from __future__ import annotations

import polars as pl


def rank_tsa_candidates(
    matrix: pl.DataFrame,
    target_tissue: str,
    top_n: int = 100,
    min_target_log2: float = 6.0,
    min_specificity_log2: float = 2.0,
) -> pl.DataFrame:
    if target_tissue not in matrix.columns:
        raise ValueError(f"Target tissue '{target_tissue}' not found in matrix.")

    tissue_cols: list[str] = [
        col for col in matrix.columns
        if col != "gene_symbol"
    ]

    other_tissues: list[str] = [
        col for col in tissue_cols
        if col != target_tissue
    ]

    if not other_tissues:
        raise ValueError("No other tissues available for comparison.")

    ranked = (
        matrix
        .filter(pl.col("gene_symbol").is_not_null())
        .with_columns(
            pl.col(target_tissue).alias("target_log2"),
            pl.max_horizontal([pl.col(col) for col in other_tissues]).alias("other_max_log2"),
            pl.mean_horizontal([pl.col(col) for col in other_tissues]).alias("other_mean_log2"),
        )
        .with_columns(
            (pl.col("target_log2") - pl.col("other_max_log2")).alias("specificity_log2"),
            (pl.col("target_log2") - pl.col("other_mean_log2")).alias("enrichment_log2"),
        )
        .with_columns(
            (
                2.0 * pl.col("specificity_log2")
                + pl.col("enrichment_log2")
                + 0.25 * pl.col("target_log2")
            ).alias("tsa_score")
        )
        .filter(pl.col("target_log2") >= min_target_log2)
        .filter(pl.col("specificity_log2") >= min_specificity_log2)
        .sort("tsa_score", descending=True)
        .head(top_n)
    )

    return ranked

def drop_poorest_duplicates(matrix: pl.DataFrame) -> pl.DataFrame:
    if "gene_symbol" not in matrix.columns:
        raise ValueError("Missing required column: gene_symbol")

    if "tsa_score" not in matrix.columns:
        raise ValueError("Missing required column: tsa_score")

    return (
        matrix
        .filter(pl.col("gene_symbol").is_not_null())
        .sort("tsa_score", descending=True)
        .group_by("gene_symbol", maintain_order=True)
        .first()
    )