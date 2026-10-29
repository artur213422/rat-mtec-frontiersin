import polars as pl


PICKED_COLUMNS: list[str] = [
    "gene_symbol",
    "chromosome",
    "strand",
    "start",
    "stop",
    "log2_fold_change",
    "p_val",
    "p_adj",
    "ortholog_id",
    "t_statistic",
    "has_clean_ortholog"
]

def cluster_genes_by_distance(dataset: pl.DataFrame, dataset_id: str, max_gene_gap_kb: int) -> pl.DataFrame:
    if dataset.is_empty():
        raise ValueError("Empty dataset")

    if dataset.select(pl.col("start").is_null().any() | pl.col("stop").is_null().any()).item():
        raise ValueError("Gene coordinates contain null values")

    max_gene_gap_bp: int = max_gene_gap_kb * 1_000

    dataset = dataset.sort(["chromosome", "start"])


    return (
        dataset
        .with_columns(
            pl.col("stop")
            .cum_max()
            .shift(1)
            .over("chromosome")
            .alias("_previous_stop")
        )
        .with_columns(
            (
                pl.col("_previous_stop").is_null()
                | (
                    pl.col("start") - pl.col("_previous_stop")
                    > max_gene_gap_bp
                )
            ).alias("_group_start")
        )
        .with_columns(
            pl.col("_group_start")
            .cast(pl.Int64)
            .cum_sum()
            .alias("_group")
        )
        .with_columns(
            pl.len()
            .over("_group")
            .alias("_group_size")
        )
        .with_columns(
            (
                (pl.col("_group") != pl.col("_group").shift(1))
                .fill_null(True)
                & (pl.col("_group_size") >= 2)
            ).alias("_cluster_start")
        )
        .with_columns(
            pl.when(pl.col("_group_size") >= 2)
            .then(
                pl.col("_cluster_start")
                .cast(pl.Int64)
                .cum_sum()
            )
            .otherwise(None)
            .alias("local_cluster_id"),
            pl.lit(dataset_id).alias("dataset_id"),
        )
        .select(["dataset_id"] + PICKED_COLUMNS + ["local_cluster_id"])
    )