from __future__ import annotations

import polars as pl


def classify_non_tsa_probes(
    dynamic_step_df: pl.DataFrame,
    *,
    probe_id_column: str = "ID_REF",
    metadata_columns: list[str] | None = None,
    tsa_column: str = "is_tissue_restricted",
) -> pl.DataFrame:
    if metadata_columns is None:
        metadata_columns = []

    if probe_id_column in metadata_columns:
        raise ValueError(f"{probe_id_column!r} must not also be included in metadata_columns")

    required_columns: set[str] = {
        probe_id_column,
        tsa_column,
        *metadata_columns,
    }

    missing_columns: set[str] = required_columns - set(dynamic_step_df.columns)

    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    if dynamic_step_df.schema[tsa_column] != pl.Boolean:
        raise TypeError(f"{tsa_column!r} must be Boolean, but its dtype is {dynamic_step_df.schema[tsa_column]}")

    tsa_null_count: int = dynamic_step_df.get_column(tsa_column).null_count()


    if tsa_null_count > 0:
        raise ValueError(f"{tsa_column!r} contains {tsa_null_count} null values")

    return (
        dynamic_step_df
        .select(
            probe_id_column,
            *metadata_columns,
            pl.col(tsa_column).alias("is_tsa"),
            (~pl.col(tsa_column)).alias("is_non_tsa"),
        )
    )