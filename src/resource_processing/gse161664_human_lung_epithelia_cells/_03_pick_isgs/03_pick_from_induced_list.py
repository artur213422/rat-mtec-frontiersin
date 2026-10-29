from pathlib import Path

import polars as pl

from src.common.feratures.multi_contditional_limma_voom.induced_specific import pick_induced_specific
from src.common.paths import find_file_in_project, resolve_dir_from_root

################################################################################
# PREPARE DATA
################################################################################

candidate_file: Path = find_file_in_project("gse161664_ifn_specific_induced_list.csv")

candidate_df: pl.DataFrame = pl.read_csv(
    candidate_file,
    infer_schema_length=None,
)

picked_isgs_file: Path = resolve_dir_from_root("resources/interferon/gse161664_lung/standardized/picked_isgs") / "gse161664_ifn_specific_induced_isgs.csv"


################################################################################
# SETTINGS
################################################################################

# Required induction above untreated control:
induction_log2_fc: float = 1.0

# Required preference over the mean response to the other IFNs:
specific_log2_fc: float = 0.5

# Applied separately to both limma contrasts:
fdr: float = 0.05


################################################################################
# SELECT INDUCED-SPECIFIC GENES
################################################################################

picked_isgs_df = pick_induced_specific(
    candidates=candidate_df,
    target_conditions=("ifna", "ifnb", "ifng", "ifnl"),
    induction_log2_fc=1.0,
    specific_log2_fc=0.5,
    fdr=0.05,
)

################################################################################
# SAVE
################################################################################

picked_isgs_df.write_csv(picked_isgs_file)

print(f"Saved: {picked_isgs_file}")