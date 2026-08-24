from pathlib import Path

import polars as pl

from src.common.paths import resolve_dir_from_root
from src.exploration.mtec_high.tsa_expression_pattern.tsa_singed_agreement_analysis.types.api_for_standardized_data import \
    AIRE_STATUS_TYPE, AIRE_GLOBAL_STATUS_TYPE

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

intermediate_dir: Path = resolve_dir_from_root("data/mtec_high/analysis/tsa_expression_pattern/tsa_singed_agreement_analysis/tsa_list/intermediate")

source_file: Path = intermediate_dir / "cleaned_matched_tsa.csv"

standardized_file: Path = resolve_dir_from_root("data/mtec_high/analysis/tsa_expression_pattern/tsa_singed_agreement_analysis/tsa_list/standardized") / "matched_tsa_list.csv"

tsa_df: pl.DataFrame = pl.read_csv(source_file)

##########################################################################################
#                                       SETTINGS
##########################################################################################

threshold_log2_foldchange: float = 0.5
threshold_p_value: float = 0.05

use_p_value: bool = True
use_raw_p_val: bool = False


##########################################################################################
#                                       MARK AIRE TSAs
##########################################################################################

p_value_col: str = "p_val" if use_raw_p_val else "p_adj"

significant: pl.Expr = (
    pl.col(p_value_col) <= threshold_p_value
    if use_p_value
    else pl.lit(True)
)

# Dataset-specific status.
tsa_df = tsa_df.with_columns(
    pl.when(
        significant
        & (pl.col("log2_fold_change") <= -threshold_log2_foldchange)
    )
    .then(pl.lit("downregulated"))
    .when(
        significant
        & (pl.col("log2_fold_change") >= threshold_log2_foldchange)
    )
    .then(pl.lit("upregulated"))
    .otherwise(pl.lit("not_changed"))
    .cast(AIRE_STATUS_TYPE)
    .alias("aire_status_dataset")
)

# Global status across datasets.
tsa_df = tsa_df.with_columns(
    pl.when(
        (pl.col("aire_status_dataset") != "not_changed")
        .any()
        .over(["tissue_group", "gene_symbol"])
    )
    .then(pl.lit("changed"))
    .otherwise(pl.lit("unchanged"))
    .cast(AIRE_GLOBAL_STATUS_TYPE)
    .alias("aire_status_global")
)
##########################################################################################
#                                       FINAL CLEANING
##########################################################################################

tsa_df = tsa_df.drop("found_in_dataset")

tsa_df = tsa_df.rename({
    "rank": "tsa_rank",
})

##########################################################################################
#                                           SAVE
##########################################################################################

tsa_df.write_csv(standardized_file)