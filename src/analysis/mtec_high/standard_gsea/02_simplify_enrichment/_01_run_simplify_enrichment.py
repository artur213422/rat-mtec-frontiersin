from __future__ import annotations

from pathlib import Path

from src.analysis.mtec_high.standard_gsea.helpers.simplify_enrichment_helpers import run_single_analyses, \
    run_comparison_analysis
from src.analysis.mtec_high.standard_gsea.helpers.apear_helpers import (
    get_standard_gsea_files,
    write_standard_gsea_direction,
)
from src.common.feratures.simplify_enrichment.simplify_enrichment import  GoOntology
from src.common.paths import PROJECT_ROOT


STANDARD_GSEA_FOLDER: Path = PROJECT_ROOT / "data/mtec_high/derived_data/standard_gsea/both_sided_gsea"
SIMPLIFY_ENRICHMENT_FOLDER: Path = PROJECT_ROOT / "data/mtec_high/derived_data/standard_gsea/simplify_enrichment"
TEMP_DIRECTIONAL_FOLDER: Path = PROJECT_ROOT / "data/temp/simplify_enrichment/standard_gsea_directional"


DIRECTION_TO_GSEA_MODE: dict[str, str] = {
    "upregulated": "pos",
    "downregulated": "neg",
}

GO_FILE_STEMS: dict[GoOntology, str] = {
    GoOntology.BP: "go_biological_process_gsea",
    GoOntology.MF: "go_molecular_function_gsea",
    GoOntology.CC: "go_cellular_component_gsea",
}


def main() -> None:

    standard_gsea_files: list[Path] = get_standard_gsea_files(STANDARD_GSEA_FOLDER)

    for direction, gsea_mode in DIRECTION_TO_GSEA_MODE.items():
        output_folder: Path = SIMPLIFY_ENRICHMENT_FOLDER / direction
        output_folder.mkdir(parents=True, exist_ok=True)

        print(f"Running simplifyEnrichment for standard GSEA {direction} (NES sign)")

        directional_files: dict[str, Path] = {
            path.stem.removesuffix(f"_{gsea_mode}"): path
            for path in write_standard_gsea_direction(
                gsea_files=standard_gsea_files,
                direction=direction,
                output_folder=TEMP_DIRECTIONAL_FOLDER,
            )
        }

        for ontology, file_stem in GO_FILE_STEMS.items():
            gsea_file: Path = directional_files[file_stem]

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


if __name__ == "__main__":
    main()
