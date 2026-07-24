from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Any


def run_rscript(
    rscript_exe: Path | str,
    r_script: Path,
    config_path: Path,
) -> None:
    command: list[str] = [
        str(rscript_exe),
        "--vanilla",
        str(r_script),
        str(config_path),
    ]

    subprocess.run(command, check=True)


def write_json_config(path: Path, config: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(config, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


project_dir = Path(r"G:\10. ÜLIÕPILASED\Artur\9. Review\analysis")

cel_files = [
    project_dir / "data/raw/cel_files/rat/Het_P20-3_mTEChi.CEL",
    project_dir / "data/raw/cel_files/rat/Het_P21-8_mTEChi.CEL",
]

config: dict[str, Any] = {
    "cel_files": [str(path) for path in cel_files],
    "sample_names": [path.stem for path in cel_files],
    "output_file": str(project_dir / "data/processed/rat_clariom_rma_matrix.txt"),
    "pd_package": "pd.clariom.s.rat.ht",
    "quote_output": True,
}

config_path = project_dir / "data/processed/rat_clariom_oligo_config.json"
r_script = project_dir / "scripts/cel_to_txt_oligo.R"

write_json_config(config_path, config)

# Use "Rscript" if it is in PATH.
# Otherwise use full Windows path, e.g.
# rscript_exe = Path(r"C:\Program Files\R\R-4.6.0\bin\Rscript.exe")
rscript_exe = "Rscript"

run_rscript(
    rscript_exe=rscript_exe,
    r_script=r_script,
    config_path=config_path,
)