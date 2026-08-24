from pathlib import Path

import polars as pl

from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.paths import try_find_file_in_project, resolve_dir_from_root

##########################################################################################
#                                       PREPARE DATA
##########################################################################################
matched_tsa_file: Path = resolve_dir_from_root("data/mtec_high/analysis/tsa_expression_pattern/tsa_singed_agreement_analysis/tsa_list/intermediate") / "matched_tsa.csv"
defined_tsa_df: pl.DataFrame = pl.read_csv(matched_tsa_file)


out_file: Path = resolve_dir_from_root("data/mtec_high/analysis/tsa_expression_pattern/tsa_singed_agreement_analysis/tsa_list/intermediate") / "matched_tsa_selected.csv"



##########################################################################################
#                                       SETTINGS
##########################################################################################

keep_cols: list[str] = [
    "dataset",
    "tissue_group",
    "has_clean_ortholog",
    "rank",
    "gene_symbol",
    "ortholog_id",
    "found_in_dataset",
    "log2_fold_change",
    "p_adj",
    "p_val",
    "t_statistic",
    "average_log2_expression",
]

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

defined_tsa_df = defined_tsa_df.select(keep_cols)

##########################################################################################
#                                           SAVE
##########################################################################################

defined_tsa_df.write_csv(out_file)



