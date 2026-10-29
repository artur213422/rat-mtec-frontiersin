from pathlib import Path

from src.analysis.mtec_high.standard_gsea.helpers.apear_helpers import validate_gsea_files
from src.analysis.mtec_high.standard_gsea.helpers.rrvgo_helpers import run_rrvgo_direction
from src.common.paths import PROJECT_ROOT

##########################################################################################
# PREPARE DATA
##########################################################################################


output_down: Path = PROJECT_ROOT / "data/mtec_high/derived_data/gsea_tsa_excluded/rrvgo/downregulated"
output_up: Path = PROJECT_ROOT / "data/mtec_high/derived_data/gsea_tsa_excluded/rrvgo/upregulated"

gsea_files_down: Path = PROJECT_ROOT / "data/mtec_high/derived_data/gsea_tsa_excluded/gsea/neg"
gsea_files_up: Path = PROJECT_ROOT / "data/mtec_high/derived_data/gsea_tsa_excluded/gsea/pos"

##########################################################################################
# SETTINGS
##########################################################################################

padj_cutoff: float = 0.05

##########################################################################################
# RUN RRVGO
##########################################################################################

downregulated_gsea_files: list[Path] = sorted(gsea_files_down.glob("*_neg.csv"))
upregulated_gsea_files: list[Path] = sorted(gsea_files_up.glob("*_pos.csv"))

validate_gsea_files(
    downregulated_gsea_files,
    gsea_mode="neg",
    go_pathway_folder=gsea_files_down,
)

validate_gsea_files(
    upregulated_gsea_files,
    gsea_mode="pos",
    go_pathway_folder=gsea_files_up,
)

run_rrvgo_direction(
    direction="downregulated",
    gsea_mode="neg",
    gsea_files=downregulated_gsea_files,
    output_folder=output_down,
    padj_cutoff=padj_cutoff,
)

run_rrvgo_direction(
    direction="upregulated",
    gsea_mode="pos",
    gsea_files=upregulated_gsea_files,
    output_folder=output_up,
    padj_cutoff=padj_cutoff,
)

print("Done")
