from __future__ import annotations

import polars as pl


def calculate_percentile(data: pl.DataFrame) -> pl.DataFrame:
    required_columns: set[str] = {
        "ctrl_log2_mean",
        "log2_fold_change",
    }

    missing_columns: set[str] = required_columns - set(data.columns)

    if missing_columns:
        raise ValueError(
            f"Missing required columns: {sorted(missing_columns)}"
        )

    if data.is_empty():
        raise ValueError(
            "Cannot calculate percentiles for an empty DataFrame."
        )

    n_valid: int = data.select(
        pl.col("ctrl_log2_mean")
        .is_not_null()
        .and_(pl.col("ctrl_log2_mean").is_finite())
        .sum()
    ).item()

    if n_valid == 0:
        return data.with_columns(
            pl.lit(None, dtype=pl.Float64).alias(
                "ctrl_expression_percentile"
            )
        )

    return data.with_columns(
        pl.when(
            pl.col("ctrl_log2_mean").is_not_null()
            & pl.col("ctrl_log2_mean").is_finite()
        )
        .then(
            100.0
            * (
                pl.col("ctrl_log2_mean").rank(
                    method="average"
                )
                - 0.5
            )
            / n_valid
        )
        .otherwise(None)
        .alias("ctrl_expression_percentile")
    )