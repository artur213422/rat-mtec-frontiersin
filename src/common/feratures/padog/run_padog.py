from dataclasses import dataclass
from pathlib import Path

from config.paths import RSCRIPT_EXE
from src.common.paths import find_file_in_project
from src.data_preprocessing.mtec_high.internal.R.run_R import run_rscript

@dataclass
class PadogInput:
    expression_matrix_file: Path
    gmt_file: Path
    output_file: Path
    control_group: tuple[str, ...]
    ko_group: tuple[str, ...]


    def to_config(self) -> dict[str, object]:
        return {
            "expression_matrix_file": str(self.expression_matrix_file),
            "gmt_file": str(self.gmt_file),
            "result_file": str(self.output_file),
            "control_group": list(self.control_group),
            "ko_group": list(self.ko_group),
        }

def run_padog(limma_input: PadogInput) -> None:
    r_script: Path = find_file_in_project("padog.R")

    run_rscript(
        rscript_exe=RSCRIPT_EXE,
        r_script=r_script,
        config=limma_input.to_config()
    )
