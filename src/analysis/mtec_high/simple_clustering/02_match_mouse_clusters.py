from pathlib import Path

import polars as pl

from src.analysis.mtec_high.simple_clustering.config import REFERENCE_DATASET
from src.analysis.mtec_high.simple_clustering.helpers.match_clusters.clean_orthologs import clean_orthologs
from src.analysis.mtec_high.simple_clustering.helpers.match_clusters.generate_mice_clustering_matching_summary import generate_mice_clustering_matching_summary
from src.analysis.mtec_high.simple_clustering.helpers.match_clusters.match_mice_clusters import match_mice_clusters
from src.common.dataset_catalog.api.get_datasets import read_file_with_ortholog_id
from src.common.paths import PROJECT_ROOT

# ----------------------------------------------------------------------------------------------------
# Data Preparation
# ----------------------------------------------------------------------------------------------------

cluster_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/genomics/simple_clustering/source/clusters.csv"
ortholog_cluster_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/genomics/simple_clustering/source/ortholog_clusters.csv"

cluster_df: pl.DataFrame = read_file_with_ortholog_id(cluster_file).filter(pl.col("dataset_id") != "rat_mtec_high")
ortholog_cluster_df: pl.DataFrame = clean_orthologs(read_file_with_ortholog_id(ortholog_cluster_file)).filter(pl.col("dataset_id") != "rat_mtec_high")

matched_cluster_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/genomics/simple_clustering/intermediate/matched_mice_clusters.csv"
matched_ortholog_cluster_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/genomics/simple_clustering/intermediate/matched_mice_ortholog_clusters.csv"

summary_file: Path = matched_cluster_file.parent / "mouse_matching_summary.txt"

# ----------------------------------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------------------------------

reference_dataset: str = REFERENCE_DATASET

# ----------------------------------------------------------------------------------------------------
# Clustering
# ----------------------------------------------------------------------------------------------------
cluster_result: pl.DataFrame = match_mice_clusters(cluster_df, reference_dataset)
ortholog_cluster_result: pl.DataFrame = match_mice_clusters(ortholog_cluster_df, reference_dataset)

# Keep unmatched genes for auditing; step 03 selects shared orthologs.

# ----------------------------------------------------------------------------------------------------
# Save
# ----------------------------------------------------------------------------------------------------


matched_cluster_file.parent.mkdir(parents=True, exist_ok=True)
cluster_result.write_csv(matched_cluster_file)
ortholog_cluster_result.write_csv(matched_ortholog_cluster_file)

summary: str = generate_mice_clustering_matching_summary(cluster_result, ortholog_cluster_result)

summary_file.write_text(summary, encoding="utf-8")