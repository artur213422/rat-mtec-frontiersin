import polars as pl

from dataclasses import dataclass
from pathlib import Path
from typing import Literal
from config.paths import RSCRIPT_EXE
from src.common.paths import find_file_in_project
from src.data_preprocessing.mtec_high.internal.R.run_R import run_rscript

@dataclass(frozen=True, slots=True)
class RrvgoInput:
    gsea_file: Path
    reference_species: str
    ontology: Literal["BP", "MF", "CC"]
    padj_cutoff: float
    output_folder: Path

    def to_config(self) -> dict[str, object]:
        return {
            "gsea_file": str(self.gsea_file),
            "reference_species": self.reference_species,
            "ontology": self.ontology,
            "padj_cutoff": self.padj_cutoff,
            "output_folder": str(self.output_folder),
        }


def run_rrvgo(rrvgo_input: RrvgoInput) -> None:
    r_script: Path = find_file_in_project("run_rrvgo.R")

    run_rscript(
        rscript_exe=RSCRIPT_EXE,
        r_script=r_script,
        config=rrvgo_input.to_config(),
    )
