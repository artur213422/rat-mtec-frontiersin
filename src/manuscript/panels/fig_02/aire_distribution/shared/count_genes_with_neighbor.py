import polars as pl


def count_genes_with_neighbor(tec_df: pl.DataFrame, *,max_distance_bp: int,) -> int:
    required_columns: set[str] = {"chromosome", "start", "stop"}

    missing_columns: set[str] = required_columns - set(tec_df.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    if tec_df.height < 2:
        raise ValueError("At least two genes are required.")

    genomic_df: pl.DataFrame = (
        tec_df
        .select(
            "chromosome",
            "start",
            "stop",
        )
        .drop_nulls()
        .sort(
            "chromosome",
            "start",
            "stop",
        )
    )

    invalid_intervals: int = genomic_df.filter(pl.col("stop") < pl.col("start")).height

    if invalid_intervals > 0:
        raise ValueError(f"Found {invalid_intervals} genes where stop < start.")

    distance_df: pl.DataFrame = (
        genomic_df
        .with_columns(
            pl.col("stop")
            .shift(1)
            .cum_max()
            .over("chromosome")
            .alias("left_max_stop"),

            pl.col("start")
            .shift(-1)
            .over("chromosome")
            .alias("right_start"),
        )
        .with_columns(
            (
                    pl.col("start") - pl.col("left_max_stop")
             )
            .clip(lower_bound=0)
            .alias("left_distance_bp"),
            (
                pl.col("right_start") - pl.col("stop")
            )
            .clip(lower_bound=0)
            .alias("right_distance_bp"),
        )
        .with_columns(
            (
                pl.col("left_distance_bp").le(max_distance_bp) | pl.col("right_distance_bp").le(max_distance_bp)
            )
            .fill_null(False)
            .alias("has_neighbor")
        )
    )

    return distance_df.filter(pl.col("has_neighbor")).height