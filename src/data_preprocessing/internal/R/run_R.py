from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any


def run_rscript(
    rscript_exe: Path,
    r_script: Path,
    config: dict[str, Any],
) -> None:
    print(f"Starting running: {r_script}")

    config_json: str = json.dumps(
        config,
        ensure_ascii=False,
        separators=(",", ":"),
    )

    command: list[str] = [
        str(rscript_exe),
        "--vanilla",
        str(r_script),
        config_json,
    ]

    subprocess.run(
        command,
        check=True,
        stdout=sys.stdout,
        stderr=sys.stderr,
    )