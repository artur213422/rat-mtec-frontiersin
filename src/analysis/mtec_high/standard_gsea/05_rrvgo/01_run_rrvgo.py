from pathlib import Path

from src.analysis.mtec_high.standard_gsea.helpers.apear_helpers import (
    get_standard_gsea_files,
    write_standard_gsea_direction,
)
from src.analysis.mtec_high.standard_gsea.helpers.rrvgo_helpers import run_rrvgo_direction

from src.common.paths import PROJECT_ROOT

##########################################################################################
# PREPARE DATA
##########################################################################################


output_up: Path = PROJECT_ROOT / "data/mtec_high/derived_data/standard_gsea/rrvgo/upregulated"
output_down: Path = PROJECT_ROOT / "data/mtec_high/derived_data/standard_gsea/rrvgo/downregulated"

gsea_folder: Path = PROJECT_ROOT / "data/mtec_high/derived_data/standard_gsea/both_sided_gsea"
directional_temp_folder: Path = PROJECT_ROOT / "data/temp/gsea/rrvgo/standard_gsea_directional"

##########################################################################################
# SETTINGS
##########################################################################################

padj_cutoff: float = 0.05

##########################################################################################
# RUN RRVGO
##########################################################################################

standard_gsea_files: list[Path] = get_standard_gsea_files(gsea_folder)
upregulated_gsea_files: list[Path] = write_standard_gsea_direction(
    gsea_files=standard_gsea_files,
    direction="upregulated",
    output_folder=directional_temp_folder,
)
downregulated_gsea_files: list[Path] = write_standard_gsea_direction(
    gsea_files=standard_gsea_files,
    direction="downregulated",
    output_folder=directional_temp_folder,
)

run_rrvgo_direction(
    direction="upregulated",
    gsea_mode="pos",
    gsea_files=upregulated_gsea_files,
    output_folder=output_up,
    padj_cutoff=padj_cutoff,
)
run_rrvgo_direction(
    direction="downregulated",
    gsea_mode="neg",
    gsea_files=downregulated_gsea_files,
    output_folder=output_down,
    padj_cutoff=padj_cutoff,
)

print("Done")
