import polars as pl

from src.analysis.mtec_high.genomics.curated_loci._01_build_physical_clusters.helpers.assign_clusters import assign_clusters
from src.analysis.mtec_high.genomics.curated_loci.utils.configs.paths import DATASET_WITH_CLUSTERS
from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata


datasets: list[DatasetMetadata] = get_mtec_high()

dataset_frames: list[pl.DataFrame] = []

for dataset in datasets:
    tec_df: pl.DataFrame = dataset.get_standardized_df()

    tec_with_clusters_df: pl.DataFrame = assign_clusters(tec_df).with_columns(pl.lit(dataset.dataset_id).alias("dataset_id"))

    dataset_frames.append(tec_with_clusters_df)

if not dataset_frames:
    raise ValueError("No datasets found.")

datasets_with_clusters_df: pl.DataFrame = pl.concat(dataset_frames, how="vertical_relaxed",)

datasets_with_clusters_df.write_csv(DATASET_WITH_CLUSTERS)

print(f"Saved {datasets_with_clusters_df.height} rows from {len(dataset_frames)} datasets to {DATASET_WITH_CLUSTERS}")