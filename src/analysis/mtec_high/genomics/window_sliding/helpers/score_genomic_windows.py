import polars as pl

from src.analysis.mtec_high.genomics.window_sliding.columns import Windows, ScoredWindows, AireDirection


def score_genomic_windows(windows_df: pl.DataFrame) -> pl.DataFrame:
    return (
        windows_df
        .with_columns(
            (
                pl.col(Windows.T_STATISTICS).list.sum()
                / pl.col(Windows.GENE_COUNT).sqrt()
            ).alias(ScoredWindows.SCORE),
        )
        .with_columns(
            pl.when(pl.col(ScoredWindows.SCORE) >= 0)
            .then(pl.lit(AireDirection.POSITIVE.value))
            .otherwise(pl.lit(AireDirection.NEGATIVE.value))
            .cast(pl.Enum(AireDirection))
            .alias(ScoredWindows.AIRE_DIRECTION),

            # Passenger handling not implemented yet.
            pl.lit([], dtype=pl.List(pl.String))
            .alias(ScoredWindows.PASSENGER_GENES),

            pl.lit(0, dtype=pl.Int64)
            .alias(ScoredWindows.PASSENGER_COUNT),

            # Initially every gene contributes to the score.
            pl.col(Windows.GENES)
            .alias(ScoredWindows.CORE_GENES),
        )
    )