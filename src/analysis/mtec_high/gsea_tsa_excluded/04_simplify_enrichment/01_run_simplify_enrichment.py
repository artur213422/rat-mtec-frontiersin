from __future__ import annotations

from pathlib import Path

import polars as pl

from src.analysis.mtec_high.standard_gsea.helpers.simplify_enrichment_helpers import run_comparison_analysis, \
    run_single_analyses
from src.common.feratures.simplify_enrichment.simplify_enrichment import GoOntology

from src.common.paths import PROJECT_ROOT


GSEA_FOLDER: Path = PROJECT_ROOT / "data/mtec_high/derived_data/gsea_tsa_excluded/gsea"
SIMPLIFY_ENRICHMENT_FOLDER: Path = PROJECT_ROOT / "data/mtec_high/derived_data/gsea_tsa_excluded/simplify_enrichment"

TEMP_FOLDER: Path = PROJECT_ROOT / "data/temp/simplify_enrichment"

# fgsea scoreType values map directly to biological direction here.
DIRECTION_TO_GSEA_MODE: dict[str, str] = {
    "upregulated": "pos",
    "down_regulated": "neg",
}

GO_FILE_STEMS: dict[GoOntology, str] = {
    GoOntology.BP: "go_biological_process_gsea",
    GoOntology.MF: "go_molecular_function_gsea",
    GoOntology.CC: "go_cellular_component_gsea",
}



for direction, gsea_mode in DIRECTION_TO_GSEA_MODE.items():
    output_folder: Path = SIMPLIFY_ENRICHMENT_FOLDER / direction
    output_folder.mkdir(parents=True, exist_ok=True)

    (output_folder / "single").mkdir(parents=True, exist_ok=True)
    (output_folder / "comparison").mkdir(parents=True, exist_ok=True)

    print(f"Running simplifyEnrichment for {direction} ({gsea_mode})")

    for ontology, file_stem in GO_FILE_STEMS.items():
        gsea_file: Path = GSEA_FOLDER / gsea_mode / f"{file_stem}_{gsea_mode}.csv"

        if not gsea_file.is_file():
            raise FileNotFoundError(f"One-sided GSEA input was not found: {gsea_file}. Run 03_run_one_sided_gsea.py first.")

        run_single_analyses(
            direction=direction,
            ontology=ontology,
            gsea_file=gsea_file,
            output_folder=output_folder,
        )


        run_comparison_analysis(
            ontology=ontology,
            gsea_file=gsea_file,
            output_folder=output_folder,
        )

    print("Done")
