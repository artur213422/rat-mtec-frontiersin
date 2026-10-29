"""Run all sliding-window genomics steps in dependency order."""

from pathlib import Path
import subprocess
import sys

from src.common.file_io import delete_all_files

PIPELINE_FOLDER = Path(__file__).resolve().parent
PROJECT_ROOT = PIPELINE_FOLDER.parents[4]
STEPS = (
    "01_build_windows.py",
    "02_score_windows.py",
    "03_select_non_overlapping_clusters.py",
    "04_drop_poor_clusters.py",
)


def main() -> None:

    delete_all_files(PROJECT_ROOT / "data/mtec_high/derived_data/genomics/window_sliding/source")
    delete_all_files(PROJECT_ROOT / "data/mtec_high/derived_data/genomics/window_sliding/standardized")

    for index, step in enumerate(STEPS, start=1):
        script = PIPELINE_FOLDER / "steps" / step
        module = ".".join(script.relative_to(PROJECT_ROOT).with_suffix("").parts)
        print(f"[{index}/{len(STEPS)}] Running: {step}", flush=True)
        subprocess.run(
            [sys.executable, "-u", "-m", module],
            cwd=PROJECT_ROOT,
            check=True,
        )
    print("Sliding-window genomics pipeline complete.", flush=True)


if __name__ == "__main__":
    main()