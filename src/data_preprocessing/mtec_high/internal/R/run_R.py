from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

def run_rscript(rscript_exe: Path, r_script: Path, config: dict[str, Any],) -> None:
    print(f"Starting running: {r_script}")

    with tempfile.NamedTemporaryFile(
            mode="w",
            suffix=".json",
            encoding="utf-8",
            delete=False,
    ) as config_file:
        json.dump(config, config_file, ensure_ascii=False, separators=(",", ":"))
        config_path = Path(config_file.name)

    try:
        subprocess.run(
            [
                str(rscript_exe),
                "--vanilla",
                str(r_script),
                str(config_path),
            ],
            check=True,
            stdout=sys.stdout,
            stderr=sys.stderr,
        )
    finally:
        config_path.unlink(missing_ok=True)