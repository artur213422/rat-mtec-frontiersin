from dataclasses import dataclass
from pathlib import Path

from config.paths import RSCRIPT_EXE
from src.common.paths import find_file_in_project
from src.data_preprocessing.mtec_high.internal.R.run_R import run_rscript


@dataclass(frozen=True, slots=True)
class APearInput:
    gsea_file: Path
    output_svg: Path
    output_clusters: Path

    def to_config(self) -> dict[str, object]:
        return {
            "gsea_file": str(self.gsea_file),
            "output_svg": str(self.output_svg),
            "output_clusters": str(self.output_clusters),
        }


def run_apear(apear_input: APearInput) -> None:
    r_script: Path = find_file_in_project("run_aPEAR.R")

    run_rscript(
        rscript_exe=RSCRIPT_EXE,
        r_script=r_script,
        config=apear_input.to_config(),
    )