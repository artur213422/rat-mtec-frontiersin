from __future__ import annotations

from pathlib import Path
from src.analysis.mtec_high.standard_gsea.helpers.apear_helpers import validate_gsea_files, run_direction
from src.common.paths import PROJECT_ROOT

##########################################################################################
# PREPARE DATA
##########################################################################################

GSEA_FOLDER: Path = PROJECT_ROOT / "data/mtec_high/derived_data/gsea_tsa_excluded/gsea"
POS_GSEA_FOLDER: Path = GSEA_FOLDER / "pos"
NEG_GSEA_FOLDER: Path = GSEA_FOLDER / "neg"
APEAR_OUTPUT_FOLDER: Path = PROJECT_ROOT / "data/mtec_high/derived_data/gsea_tsa_excluded/apear"
TEMP_FOLDER: Path = PROJECT_ROOT / "data/temp/gsea/apear"

GO_FILE_STEM_TO_ABBREVIATION: dict[str, str] = {
    "go_biological_process_gsea": "BP",
    "go_cellular_component_gsea": "CC",
    "go_molecular_function_gsea": "MF",
}

# fgsea scoreType values map directly to biological direction here.
upregulated_gsea_files: list[Path] = sorted(POS_GSEA_FOLDER.glob("*_pos.csv"))
downregulated_gsea_files: list[Path] = sorted(NEG_GSEA_FOLDER.glob("*_neg.csv"))



validate_gsea_files(upregulated_gsea_files, gsea_mode="pos", go_pathway_folder=POS_GSEA_FOLDER)
validate_gsea_files(downregulated_gsea_files, gsea_mode="neg", go_pathway_folder=NEG_GSEA_FOLDER)

run_direction(
    direction="upregulated",
    gsea_mode="pos",
    gsea_files=upregulated_gsea_files,
    apear_input_folder=APEAR_OUTPUT_FOLDER
)

run_direction(
    direction="down_regulated",
    gsea_mode="neg",
    gsea_files=downregulated_gsea_files,
    apear_input_folder=APEAR_OUTPUT_FOLDER
)

print("Done")
