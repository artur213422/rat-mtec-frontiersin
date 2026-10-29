from pathlib import Path

import polars as pl

from src.common.dataset_catalog.api.get_datasets import read_file_with_ortholog_id
from src.common.paths import find_file_in_project

# -------------------------------------------------------------------------------------------------------------
# Prepare Data
# -------------------------------------------------------------------------------------------------------------

loci_file: Path = find_file_in_project("datasets_with_clusters.csv")
loci_df: pl.DataFrame = read_file_with_ortholog_id(loci_file)

cleaned_datasets: list[pl.DataFrame] = []

# -------------------------------------------------------------------------------------------------------------
# Remove Weird Chromosomes
# -------------------------------------------------------------------------------------------------------------

for tec_df in loci_df.partition_by("dataset_id", maintain_order=True):
    dataset_id: str = tec_df.get_column("dataset_id").item(0)
    n_before: int = tec_df.height

    valid_chromosome = pl.col("chromosome").str.contains(r"^chr(?:[0-9]+|X|Y)$").fill_null(False)

    removed_chromosomes: list[str] = (
        tec_df
        .filter(~valid_chromosome)
        .get_column("chromosome")
        .drop_nulls()
        .unique()
        .sort()
        .to_list()
    )

    tec_df = tec_df.filter(valid_chromosome)

    n_after: int = tec_df.height
    n_dropped: int = n_before - n_after
    dropped_percent: float = (
        (n_dropped / n_before) * 100 if n_before > 0 else 0.0
    )

    # ---------------------------------------------------------------------------------------------------------
    # Re-index clusters continuously within this dataset: 0, 1, 2, ...
    # ---------------------------------------------------------------------------------------------------------

    cluster_map_df: pl.DataFrame = (
        tec_df
        .group_by("cluster_id")
        .agg(
            pl.col("chromosome").first().alias("chromosome"),
            pl.col("start").min().alias("cluster_start"),
        )
        .with_columns(
            pl.when(pl.col("chromosome") == "chrX")
            .then(pl.lit(1_000))
            .when(pl.col("chromosome") == "chrY")
            .then(pl.lit(1_001))
            .otherwise(
                pl.col("chromosome")
                .str.extract(r"^chr(\d+)$", 1)
                .cast(pl.Int64)
            )
            .alias("_chromosome_order")
        )
        .sort(
            "_chromosome_order",
            "cluster_start",
            "cluster_id",
        )
        .with_row_index(name="new_cluster_id", offset=0)
        .select("cluster_id", "new_cluster_id")
    )

    tec_df = (
        tec_df
        .join(cluster_map_df, on="cluster_id", how="left")
        .drop("cluster_id")
        .rename({"new_cluster_id": "cluster_id"})
        .select(loci_df.columns)
        .sort("cluster_id", "start", "stop")
    )

    cleaned_datasets.append(tec_df)

    # ---------------------------------------------------------------------------------------------------------
    # Report
    # ---------------------------------------------------------------------------------------------------------

    n_clusters: int = tec_df.get_column("cluster_id").n_unique()

    print(f"{dataset_id}: dropped {n_dropped}/{n_before} genes ({dropped_percent:.2f}%). {n_after} genes remain in {n_clusters} clusters.")

    if removed_chromosomes:
        print(f"  Removed chromosomes: {', '.join(removed_chromosomes)}")

# -------------------------------------------------------------------------------------------------------------
# Save All Datasets to One File
# -------------------------------------------------------------------------------------------------------------

cleaned_loci_df: pl.DataFrame = (
    pl.concat(cleaned_datasets, how="vertical")
    if cleaned_datasets
    else loci_df
)

cleaned_loci_df.write_csv(loci_file)

print(f"Saved {cleaned_loci_df.height} rows to {loci_file}")