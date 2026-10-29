from __future__ import annotations

import polars as pl

from src.analysis.mtec_high.genomics.window_sliding.columns import ScoredWindows, Windows


def select_non_overlapping_clusters(windows_df: pl.DataFrame) -> pl.DataFrame:
    selected_chromosomes: list[pl.DataFrame] = []

    for chromosome_df in windows_df.partition_by(Windows.CHROMOSOME, maintain_order=True):
        chromosome_df = chromosome_df.with_row_index("_row_id")

        ranked_df: pl.DataFrame = (
            chromosome_df
            .with_columns(
                pl.col(ScoredWindows.SCORE)
                .abs()
                .alias("_absolute_score")
            )
            .sort("_absolute_score", descending=True)
        )

        selected_row_ids: list[int] = []
        selected_intervals: list[tuple[int, int]] = []

        for row in ranked_df.iter_rows(named=True):
            window_start: int = row[Windows.START]
            window_stop: int = row[Windows.STOP]

            overlaps_existing: bool = any(
                window_start <= selected_stop
                and window_stop >= selected_start
                for selected_start, selected_stop in selected_intervals
            )

            if overlaps_existing:
                continue

            selected_row_ids.append(row["_row_id"])
            selected_intervals.append(
                (window_start, window_stop)
            )

        selected_df: pl.DataFrame = (
            chromosome_df
            .filter(pl.col("_row_id").is_in(selected_row_ids))
            .drop("_row_id")
            .sort(Windows.START, Windows.STOP)
        )

        selected_chromosomes.append(selected_df)

    if not selected_chromosomes:
        return windows_df.head(0)

    return pl.concat(selected_chromosomes, how="vertical")