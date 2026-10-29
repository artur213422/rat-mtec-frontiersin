from enum import StrEnum

from dataclasses import dataclass
from pathlib import Path
from config.paths import RSCRIPT_EXE
from src.common.paths import find_file_in_project
from src.common.types.enums import Species
from src.data_preprocessing.mtec_high.internal.R.run_R import run_rscript

class GoOntology(StrEnum):
    BP = "BP"
    CC = "CC"
    MF = "MF"

class SimplifyEnrichmentMode(StrEnum):
    SINGLE = "single"
    COMPARISON = "comparison"

@dataclass(frozen=True, slots=True)
class SimplifyEnrichmentInput:
    gsea_file: Path
    output_svg: Path
    output_clusters: Path

    ontology: GoOntology
    reference_species: Species
    mode: SimplifyEnrichmentMode

    padj_cutoff: float = 0.05

    def __post_init__(self) -> None:
        if not 0.0 < self.padj_cutoff <= 1.0:
            raise ValueError(f"padj_cutoff must be in (0, 1], got {self.padj_cutoff}")

    def to_config(self) -> dict[str, object]:
        return {
            "gsea_file": str(self.gsea_file),
            "output_svg": str(self.output_svg),
            "output_clusters": str(self.output_clusters),
            "ontology": self.ontology.value,
            "reference_species": self.reference_species.value,
            "analysis_mode": self.mode.value,
            "padj_cutoff": self.padj_cutoff,
        }



def run_simplify_enrichment(simplify_enrichment_input: SimplifyEnrichmentInput) -> None:
    r_script: Path = find_file_in_project("simplify_enrichment.R")

    run_rscript(
        rscript_exe=RSCRIPT_EXE,
        r_script=r_script,
        config=simplify_enrichment_input.to_config(),
    )
