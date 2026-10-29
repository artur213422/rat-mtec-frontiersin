from pathlib import Path

import polars as pl

from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.paths import resolve_dir_from_root, find_file_in_project, PROJECT_ROOT

##########################################################################################
# PREPARE DATA
##########################################################################################

aps1_antigenes_file: Path = find_file_in_project("aps1_antigens_with_resolved_orthologs.csv")

aps1_antigenes_df: pl.DataFrame = pl.read_csv(aps1_antigenes_file)

out_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/tsa_expression/aps1_known_antigens/standardized/aps1_antigens.csv"

datasets: list[DatasetMetadata] = get_mtec_high()


##########################################################################################
# Handle
##########################################################################################

stats_columns = [
    "ortholog_id",
    "gene_symbol",
    "log2_fold_change",
    "ranking_statistic",
    "p_val",
    "p_adj",
    "average_log2_expression",
]

tissue_dfs = []

for dataset in datasets:
    dataset_df: pl.DataFrame = dataset.get_standardized_df()

    stats_df = (
        dataset_df
        .with_columns(pl.col("ortholog_id").cast(pl.String))
        .select(stats_columns)
    )

    tissue_df = (
        aps1_antigenes_df
        .with_columns(
            pl.col("ortholog_id").cast(pl.String),
            pl.lit(dataset.dataset_id).alias("dataset_id"),
        )
        .join(stats_df, on="ortholog_id", how="left")
    )

    tissue_dfs.append(tissue_df)

tissue_df = pl.concat(tissue_dfs)



##########################################################################################
# SAVE
##########################################################################################

tissue_df.write_csv(out_file)