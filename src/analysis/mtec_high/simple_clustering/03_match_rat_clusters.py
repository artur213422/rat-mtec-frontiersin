from pathlib import Path

import polars as pl

from src.analysis.mtec_high.simple_clustering.config import REFERENCE_DATASET
from src.analysis.mtec_high.simple_clustering.helpers.match_clusters.match_rat_clusters import match_rat_clusters, select_shared_orthologs
from src.common.dataset_catalog.api.get_datasets import read_file_with_ortholog_id
from src.common.paths import PROJECT_ROOT

# ----------------------------------------------------------------------------------------------------
# Data Preparation
# ----------------------------------------------------------------------------------------------------

base: Path = PROJECT_ROOT / "data/mtec_high/derived_data/genomics/simple_clustering"

mice: pl.DataFrame = read_file_with_ortholog_id(base / "intermediate/matched_mice_ortholog_clusters.csv")
rat: pl.DataFrame = read_file_with_ortholog_id(base / "source/ortholog_clusters.csv").filter(pl.col("dataset_id") == "rat_mtec_high")

output: Path = base / "standardized/matched_clusters.csv"

# ----------------------------------------------------------------------------------------------------
# Matching
# ----------------------------------------------------------------------------------------------------

result: pl.DataFrame = match_rat_clusters(mice, rat, REFERENCE_DATASET)
shared: pl.DataFrame = select_shared_orthologs(result, REFERENCE_DATASET)

# ----------------------------------------------------------------------------------------------------
# Save
# ----------------------------------------------------------------------------------------------------

output.parent.mkdir(parents=True, exist_ok=True)
shared.write_csv(output)

print(f"Shared ortholog pairs: {shared.filter(pl.col('dataset_id') == REFERENCE_DATASET).height}")