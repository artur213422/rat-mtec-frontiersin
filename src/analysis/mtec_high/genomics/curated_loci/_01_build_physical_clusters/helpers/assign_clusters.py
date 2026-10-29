import polars as pl

from src.analysis.mtec_high.genomics.curated_loci.utils.configs.config import CLUSTER_GAP_KB, PHYSICAL_CLUSTER_COLS


def assign_clusters(tec_df: pl.DataFrame) -> pl.DataFrame:
    required_columns: set[str] = {
        "gene_symbol",
        "chromosome",
        "start",
        "stop",
    }

    missing_columns: set[str] = required_columns - set(tec_df.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    cluster_gap_bp: int = CLUSTER_GAP_KB * 1_000

    # ---------------------------------------------------------------------------------------------------------
    # Prepare genomic data
    # ---------------------------------------------------------------------------------------------------------

    genomic_df: pl.DataFrame = (
        tec_df
        .filter(
            pl.col("chromosome").is_not_null(),
            pl.col("start").is_not_null(),
            pl.col("stop").is_not_null(),
        )
        .sort(
            "chromosome",
            "start",
            "stop",
        )
    )

    invalid_intervals: int = genomic_df.filter(pl.col("stop") < pl.col("start")).height

    if invalid_intervals > 0:
        raise ValueError(f"Found {invalid_intervals} genes where stop < start.")

    # ---------------------------------------------------------------------------------------------------------
    # Assign physical clusters
    # ---------------------------------------------------------------------------------------------------------

    clustered_df: pl.DataFrame = (
        genomic_df
        .with_columns(
            pl.col("stop")
            .cum_max()
            .shift(1)
            .over("chromosome")
            .alias("_previous_max_stop")
        )
        .with_columns(
            (
                pl.col("start") - pl.col("_previous_max_stop")
            )
            .clip(lower_bound=0)
            .alias("_previous_gene_gap_bp")
        )
        .with_columns(
            (
                pl.col("_previous_max_stop").is_null()
                | (
                    pl.col("_previous_gene_gap_bp") > cluster_gap_bp
                )
            )
            .alias("_starts_new_cluster")
        )
        .with_columns(
            (
                pl.col("_starts_new_cluster")
                .cast(pl.Int64)
                .cum_sum()
                - 1
            )
            .alias("cluster_id")
        )
    )

    # ---------------------------------------------------------------------------------------------------------
    # Calculate cluster-level properties
    # ---------------------------------------------------------------------------------------------------------

    cluster_summary_df: pl.DataFrame = (
        clustered_df
        .group_by("cluster_id", maintain_order=True)
        .agg(
            pl.col("chromosome").first().alias("chromosome"),
            pl.col("start").min().alias("cluster_start"),
            pl.col("stop").max().alias("cluster_stop"),
            pl.len().alias("cluster_gene_count"),
        )
        .sort("chromosome", "cluster_start",)
        .with_columns(
            pl.col("cluster_stop")
            .shift(1)
            .over("chromosome")
            .alias("_previous_cluster_stop"),

            pl.col("cluster_start")
            .shift(-1)
            .over("chromosome")
            .alias("_next_cluster_start"),
        )
        .with_columns(
            (
                pl.col("cluster_start")
                - pl.col("_previous_cluster_stop")
            )
            .clip(lower_bound=0)
            .alias("gap_to_previous_cluster_bp"),

            (
                pl.col("_next_cluster_start")
                - pl.col("cluster_stop")
            )
            .clip(lower_bound=0)
            .alias("gap_to_next_cluster_bp"),
        )
        .drop(
            "_previous_cluster_stop",
            "_next_cluster_start",
        )
    )

    # ---------------------------------------------------------------------------------------------------------
    # Add cluster properties to every gene
    # ---------------------------------------------------------------------------------------------------------

    result_df: pl.DataFrame = (
        clustered_df
        .join(
            cluster_summary_df.select(
                "cluster_id",
                "cluster_gene_count",
                "gap_to_previous_cluster_bp",
                "gap_to_next_cluster_bp",
            ),
            on="cluster_id",
            how="left",
        )
        .drop(
            "_previous_max_stop",
            "_previous_gene_gap_bp",
            "_starts_new_cluster",
        )
        .sort(
            "cluster_id",
            "start",
            "stop",
        )
    )

    available_selected_cols: list[str] = [
        column
        for column in PHYSICAL_CLUSTER_COLS
        if column in result_df.columns
    ]

    return result_df.select(available_selected_cols)