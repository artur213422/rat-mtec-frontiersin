from pathlib import Path

import polars as pl

from src.common.feratures.simplify_enrichment.simplify_enrichment import (
    GoOntology,
    SimplifyEnrichmentInput,
    SimplifyEnrichmentMode,
    run_simplify_enrichment,
)
from src.common.paths import PROJECT_ROOT

from src.common.types.enums import Species


TEMP_FOLDER: Path = PROJECT_ROOT / "data/temp/simplify_enrichment"


def get_reference_species(dataset_id: str) -> Species:
    return Species.RAT if "rat" in dataset_id.lower() else Species.MOUSE


def run_single_analyses(
    *,
    direction: str,
    ontology: GoOntology,
    gsea_file: Path,
    output_folder: Path,
) -> None:
    TEMP_FOLDER.mkdir(parents=True, exist_ok=True)

    gsea_df: pl.DataFrame = pl.read_csv(gsea_file)

    for (dataset_id_value,), dataset_df in gsea_df.group_by("dataset_id",maintain_order=True,):
        dataset_id: str = str(dataset_id_value)
        ontology_name: str = ontology.value.lower()
        output_prefix: str = f"{dataset_id}_{ontology_name}"
        dataset_file: Path = (
            TEMP_FOLDER / f"{dataset_id}_{ontology_name}_{direction}.csv"
        )

        dataset_df.write_csv(dataset_file)

        run_simplify_enrichment(
            SimplifyEnrichmentInput(
                gsea_file=dataset_file,
                output_svg=(
                    output_folder / "single" / f"{output_prefix}_simplified_enrichment.svg"
                ),
                output_clusters=output_folder / "single" / f"{output_prefix}_clusters.csv",
                ontology=ontology,
                reference_species=get_reference_species(dataset_id),
                mode=SimplifyEnrichmentMode.SINGLE,
            )
        )


def run_comparison_analysis(
    *,
    ontology: GoOntology,
    gsea_file: Path,
    output_folder: Path,
) -> None:
    ontology_name: str = ontology.value.lower()

    run_simplify_enrichment(
        SimplifyEnrichmentInput(
            gsea_file=gsea_file,
            output_svg=output_folder /"comparison" / f"{ontology_name}_simplified_enrichment.svg",
            output_clusters=output_folder /"comparison"  / f"{ontology_name}_clusters.csv",
            ontology=ontology,
            reference_species=Species.MOUSE,
            mode=SimplifyEnrichmentMode.COMPARISON,
        )
    )
