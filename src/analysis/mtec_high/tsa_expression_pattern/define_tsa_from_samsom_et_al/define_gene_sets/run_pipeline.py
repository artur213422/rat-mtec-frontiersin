"""Run atlas gene-set preparation in dependency order."""
from pathlib import Path
import subprocess
import sys

PIPELINE_FOLDER = Path(__file__).resolve().parent
STEPS = (
    "_01_define_groups/01_define_physiological_groups.py",
    "_02_define_tsa/01_run_dynamic_step_classification.py",
    "_02_define_tsa/02_final_gene_cleanup.py",
    "_03_define_non_tsa/01_define_non_tsa.py",
    "_03_define_non_tsa/02_final_gene_cleanup.py",
    "_04_match_with_datasets/01_match_with_datasets.py",
)


def main() -> None:
    for step in STEPS:
        script = PIPELINE_FOLDER / step
        print(f"Running: {step}", flush=True)
        subprocess.run([sys.executable, str(script)], check=True)
    print("Gene-set pipeline complete.")


if __name__ == "__main__":
    main()
