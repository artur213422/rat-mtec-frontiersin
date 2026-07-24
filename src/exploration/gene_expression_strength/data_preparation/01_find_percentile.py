import polars as pl
from pathlib import Path

from src.common.feratures.math.fold_change import calculate_fold_changes
from src.common.feratures.math.calc_percentile import calculate_percentile
from src.common.feratures.math.p_val import calculate_log2_mean
from src.common.models.dataset_metadata.model import DatasetCatalog
from src.common.models.dataset_metadata.ops import build_dataset_catalog
from src.common.paths import resolve_dir_from_root, find_file_in_project

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


output_file: Path = resolve_dir_from_root("data/analysis/gene_expression_strength/data")

datasets: DatasetCatalog = build_dataset_catalog()

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

for dataset in datasets.datasets:
    tec_df: pl.DataFrame = dataset.get_standardized_df()

    tec_df = calculate_fold_changes(tec_df)

    tec_df = calculate_percentile(tec_df)

    ##########################################################################################
    #                                       SAVE
    ##########################################################################################

    file: Path = output_file / f"{dataset.dataset_id}.csv"

    tec_df.write_csv(file)



