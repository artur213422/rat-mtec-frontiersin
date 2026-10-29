from pathlib import Path

import polars as pl

from src.common.feratures.multi_contditional_limma_voom.induced_specific import pick_induced_specific
from src.common.paths import resolve_dir_from_root, find_file_in_project

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

candidate_file: Path = find_file_in_project("gse28183_ifn_specific_induced_list.csv")
candidate_df: pl.DataFrame = pl.read_csv(
    candidate_file,
    infer_schema_length=None,
)

picked_isgs_file: Path = resolve_dir_from_root("resources/interferon/gse28183_tec/standardized/picked_isgs") / "gse28183_ifn_specific_induced_isgs.csv"

##########################################################################################
#                                       SETTINGS
##########################################################################################

induction_log2_fc: float = 1.0
specific_log2_fc: float = 0.5
fdr: float = 0.05


##########################################################################################
#                                       PROCESS DATA
##########################################################################################

picked_isgs_df = pick_induced_specific(
    candidates=candidate_df,
    target_conditions=("ifna", "ifng"),
    induction_log2_fc=1.0,
    specific_log2_fc=0.5,
    fdr=0.05,
)

##########################################################################################
#                                      SAVE
##########################################################################################
picked_isgs_df.write_csv(picked_isgs_file)