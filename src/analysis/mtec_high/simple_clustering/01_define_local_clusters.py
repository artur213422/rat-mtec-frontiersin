from pathlib import Path

import polars as pl

from src.analysis.mtec_high.simple_clustering.config import MAX_GENE_GAP_KB
from src.analysis.mtec_high.simple_clustering.helpers.define_local_clusters.define_cluster import cluster_genes_by_distance
from src.analysis.mtec_high.simple_clustering.helpers.define_local_clusters.generate_clustering_summary import generate_clustering_summary
from src.analysis.mtec_high.simple_clustering.helpers.define_local_clusters.get_primary_chromosome import get_primary_chromosome_df, get_primary_chromosome_orthologs_df
from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.paths import PROJECT_ROOT

# ----------------------------------------------------------------------------------------------------
# Data Preparation
# ----------------------------------------------------------------------------------------------------

datasets: list[DatasetMetadata] = get_mtec_high()

cluster_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/genomics/simple_clustering/source/clusters.csv"
ortholog_cluster_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/genomics/simple_clustering/source/ortholog_clusters.csv"
summary_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/genomics/simple_clustering/source/cluster_summary.txt"

# ----------------------------------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------------------------------

max_gene_gap_kb: int = MAX_GENE_GAP_KB

# ----------------------------------------------------------------------------------------------------
# Clustering
# ----------------------------------------------------------------------------------------------------

clusters: list[pl.DataFrame] =[]
ortholog_clusters: list[pl.DataFrame] =[]

for dataset in datasets:
    tec_df: pl.DataFrame = get_primary_chromosome_df(dataset)
    ortholog_tec_df: pl.DataFrame = get_primary_chromosome_orthologs_df(dataset)

    clusters_df: pl.DataFrame = cluster_genes_by_distance(tec_df, dataset.dataset_id, max_gene_gap_kb)
    ortholog_clusters_df: pl.DataFrame = cluster_genes_by_distance(ortholog_tec_df, dataset.dataset_id, max_gene_gap_kb)

    clusters.append(clusters_df)
    ortholog_clusters.append(ortholog_clusters_df)

# ----------------------------------------------------------------------------------------------------
# Save
# ----------------------------------------------------------------------------------------------------

cluster_result: pl.DataFrame = pl.concat(clusters)
ortholog_cluster_result: pl.DataFrame = pl.concat(ortholog_clusters)

cluster_file.parent.mkdir(parents=True, exist_ok=True)
cluster_result.write_csv(cluster_file)
ortholog_cluster_result.write_csv(ortholog_cluster_file)

summary: str = generate_clustering_summary(cluster_result, ortholog_cluster_result, max_gene_gap_kb)
summary_file.write_text(summary, encoding="utf-8")