import re

import polars as pl

def extract_tissue_replicates(tissue: str, matrix: pl.DataFrame) -> list[str]:
    """
       Return the two replicate columns for a tissue.

       Example:
           tissue="liver"
           -> ["liver_1", "liver_2"]
       """


    pattern = re.compile(rf"^{re.escape(tissue)}_\d+$")

    replicates = [
        col
        for col in matrix.columns
        if pattern.fullmatch(col)
    ]

    if len(replicates) != 2:
        raise ValueError(
            f"Expected exactly 2 replicate columns for '{tissue}', "
            f"found {len(replicates)}: {replicates}"
        )

    return sorted(replicates)

def log2_matrix(matrix: pl.DataFrame) -> pl.DataFrame:
    """Apply log2(x + 1) to all expression columns, leaving gene_symbol unchanged."""

    return matrix.with_columns(
        [
            (pl.col(col).cast(pl.Float64) + 1.0).log(base=2)
            for col in matrix.columns
            if col != "gene_symbol"
        ]
    )

def collapse_cols(matrix: pl.DataFrame, tissue: str, tissue_replicates: list[str]) -> pl.DataFrame:
    if len(tissue_replicates) != 2:
        raise ValueError(
            f"Expected exactly 2 replicate columns for tissue '{tissue}', "
            f"found {len(tissue_replicates)}: {tissue_replicates}"
        )

    missing_cols: list[str] = [
        col for col in tissue_replicates
        if col not in matrix.columns
    ]

    if missing_cols:
        raise ValueError(
            f"Cannot collapse tissue '{tissue}'. "
            f"These replicate columns are missing from matrix: {missing_cols}. "
            f"Available columns: {matrix.columns}"
        )

    return (
        matrix
        .with_columns(
            pl.mean_horizontal(
                [pl.col(col).cast(pl.Float64) for col in tissue_replicates]
            ).alias(tissue)
        )
        .drop(tissue_replicates)
    )