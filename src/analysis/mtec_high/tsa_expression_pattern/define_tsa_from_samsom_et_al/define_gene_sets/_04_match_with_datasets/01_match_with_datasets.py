from pathlib import Path
import sys

# Support direct execution from any working directory.
PROJECT_ROOT = Path(__file__).resolve().parents[7]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.analysis.mtec_high.tsa_expression_pattern.define_tsa_from_samsom_et_al.define_gene_sets.config.paths import (
    TSA_FILE,
    MATCHED_TSA_FILE,
    STANDARDIZED_FOLDER,
)

import polars as pl

from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata


##########################################################################################
# PREPARE DATA
##########################################################################################

atlas_file: Path = TSA_FILE

atlas_df: pl.DataFrame = pl.read_csv(atlas_file)

out_file: Path = MATCHED_TSA_FILE
out_file.parent.mkdir(parents=True, exist_ok=True)

datasets: list[DatasetMetadata] = get_mtec_high()


##########################################################################################
# Handle
##########################################################################################

stats_columns = [
    "gene_symbol",
    "log2_fold_change",
    "ranking_statistic",
    "p_val",
    "p_adj",
    "average_log2_expression",
]

tissue_dfs = []

for dataset in datasets:
    dataset_df: pl.DataFrame = dataset.get_standardized_orthologs_df(root=STANDARDIZED_FOLDER)

    stats_df = dataset_df.select(stats_columns)

    tissue_df = (
        atlas_df
        .with_columns(pl.lit(dataset.dataset_id).alias("dataset_id"))
        .join(stats_df, on="gene_symbol", how="left",)
    )

    tissue_dfs.append(tissue_df)

tissue_df = pl.concat(tissue_dfs)



##########################################################################################
# SAVE
##########################################################################################

tissue_df.write_csv(out_file)