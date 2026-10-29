from pathlib import Path

import polars as pl

from src.common.paths import resolve_dir_from_root, find_file_in_project

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

candidate_file: Path = find_file_in_project("gse161664_ifn_vs_control_list.csv")
candidate_df: pl.DataFrame = pl.read_csv(candidate_file)

picked_isgs_file: Path = resolve_dir_from_root("resources/interferon/gse161664_lung/standardized/picked_isgs") / "gse161664_ifn_vs_control_isgs.csv"

##########################################################################################
#                                       SETTINGS
##########################################################################################

log2_fold_change: float = 1
p_adj: float = 0.05

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

picked_isgs_df = (
    candidate_df
    # Pick DEGs
    .filter(
        (pl.col("log2_fold_change") >= log2_fold_change)
        & (pl.col("p_adj") <= p_adj)
        & pl.col("has_clear_ortholog")
        & pl.col("ortholog_id").is_not_null()
    )
    # Rank duplicate ortholog mappings:
    # lowest adjusted p-value, then largest absolute FC
    .with_columns(
        pl.col("log2_fold_change").abs().alias("_abs_log2_fc")
    )
    .sort(
        ["target_condition", "ortholog_id", "p_adj", "_abs_log2_fc"],
        descending=[False, False, False, True],
    )
    # One ortholog per IFN
    .unique(
        subset=["target_condition", "ortholog_id"],
        keep="first",
        maintain_order=True,
    )
    .drop("_abs_log2_fc")
)


##########################################################################################
#                                      SAVE
##########################################################################################
picked_isgs_df.write_csv(picked_isgs_file)