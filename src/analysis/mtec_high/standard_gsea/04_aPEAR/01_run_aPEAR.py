from __future__ import annotations

from pathlib import Path

from src.analysis.mtec_high.standard_gsea.helpers.apear_helpers import (
    get_standard_gsea_files,
    run_direction,
    write_standard_gsea_direction,
)
from src.common.paths import PROJECT_ROOT

##########################################################################################
# PREPARE DATA
##########################################################################################

GO_PATHWAYS_FOLDER: Path = PROJECT_ROOT / "data/mtec_high/derived_data/standard_gsea/both_sided_gsea"
APEAR_OUTPUT_FOLDER: Path = PROJECT_ROOT / "data/mtec_high/derived_data/standard_gsea/apear"
TEMP_FOLDER: Path = PROJECT_ROOT / "data/temp/gsea/apear/standard_gsea_directional"



standard_gsea_files: list[Path] = get_standard_gsea_files(GO_PATHWAYS_FOLDER)
upregulated_gsea_files: list[Path] = write_standard_gsea_direction(
    gsea_files=standard_gsea_files,
    direction="upregulated",
    output_folder=TEMP_FOLDER,
)
downregulated_gsea_files: list[Path] = write_standard_gsea_direction(
    gsea_files=standard_gsea_files,
    direction="downregulated",
    output_folder=TEMP_FOLDER,
)

run_direction(
    direction="upregulated",
    gsea_mode="pos",
    gsea_files=upregulated_gsea_files,
    apear_input_folder=APEAR_OUTPUT_FOLDER
)

run_direction(
    direction="downregulated",
    gsea_mode="neg",
    gsea_files=downregulated_gsea_files,
    apear_input_folder=APEAR_OUTPUT_FOLDER
)

print("Done")
