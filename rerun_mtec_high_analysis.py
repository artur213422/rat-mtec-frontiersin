from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from datetime import datetime


ROOT = Path(__file__).resolve().parent
BASE = ROOT / "src/analysis/mtec_high"
WORK = ROOT / "data/temp/mtec_high_master"
STATE_FILE = WORK / "state.json"

JOBS = {}


def add(name, folders, dependencies=(), guard=None):
    """Folders are processed in the supplied dependency order."""
    JOBS[name] = {
        "folders": folders.splitlines(),
        "dependencies": tuple(dependencies),
        "guard": guard,
    }


# Each folder contributes only its immediate numbered Python scripts.
# Helpers, tests and run_pipeline.py wrappers are not executed.
add("standard_reactome", "standard_gsea/01_gsea")
# Split this folder into independently runnable analyses below.
JOBS["standard_reactome"]["files"] = [
    "standard_gsea/01_gsea/01_reactome_gsea.py"
]
add("standard_go", "")
JOBS["standard_go"]["files"] = [
    "standard_gsea/01_gsea/02_go_gsea.py"
]
add("standard_one_sided", "")
JOBS["standard_one_sided"]["files"] = [
    "standard_gsea/01_gsea/03_run_one_sided_gsea.py"
]

# These current implementations consume two-sided GO results.
add(
    "standard_simplify",
    "standard_gsea/02_simplify_enrichment",
    ("standard_go",),
)
add("standard_apear", "standard_gsea/04_aPEAR", ("standard_go",))
add("standard_rrvgo", "standard_gsea/05_rrvgo", ("standard_go",))

add(
    "isg_tec",
    "isg_pathways_gsea/gse281833_human_tec",
)
add("isg_lung", "")
JOBS["isg_lung"]["files"] = [
    "isg_pathways_gsea/gse161664_human_lung_epithelia_cells/01_gsea.py"
]
add("isg_padog", "")
JOBS["isg_padog"]["files"] = [
    "isg_pathways_gsea/gse161664_human_lung_epithelia_cells/02_run_padog.py"
]

ATLAS = (
    "tsa_expression_pattern/"
    "define_tsa_from_samsom_et_al"
)
add(
    "tsa_atlas",
    "\n".join(
        f"{ATLAS}/define_gene_sets/{folder}"
        for folder in (
            "_01_define_groups",
            "_02_define_tsa",
            "_03_define_non_tsa",
            "_04_match_with_datasets",
        )
    ),
)
add(
    "tsa_tissue",
    f"{ATLAS}/tissue_gsea",
    ("tsa_atlas",),
)
add(
    "aps1",
    "tsa_expression_pattern/define_known_aps1_antigens",
)

add("chemokines", "chemokine_expression")
add("kallikreins", "kallikrein_expression")
add("mucins", "muc_expression")
add("mimetic", "mimetic_cells")
add("hollander", "hollander_thymus_cell_marker_gsea")
add("antigen_presentation", "antigen_presentation_gsea")
add("peptide_hormones", "peptide_hormone_gsea")

# Execute steps directly, avoiding the wrapper's delete-before-run behavior.
add("genomic_windows", "genomics/window_sliding/steps")
add("simple_clustering", "simple_clustering")
add(
    "curated_loci",
    "\n".join(
        (
            "genomics/curated_loci/_01_build_physical_clusters",
            "genomics/curated_loci/_02_map_curated_loci/gene_set",
            "genomics/curated_loci/_02_map_curated_loci/gene_family",
        )
    ),
)

add(
    "tsa_excluded_prepare",
    "gsea_tsa_excluded/01_prepare_tsa_excluded_datasets",
    ("tsa_atlas",),
)
add(
    "tsa_excluded_go",
    "gsea_tsa_excluded/02_go_gsea",
    ("tsa_excluded_prepare",),
)
add(
    "tsa_excluded_apear",
    "gsea_tsa_excluded/03_apear",
    ("tsa_excluded_go",),
    guard="tsa_directional",
)
add(
    "tsa_excluded_simplify",
    "gsea_tsa_excluded/04_simplify_enrichment",
    ("tsa_excluded_go",),
    guard="tsa_directional",
)
add(
    "tsa_excluded_rrvgo",
    "gsea_tsa_excluded/05_rrvgo",
    ("tsa_excluded_go",),
    guard="tsa_directional",
)


def scripts_for(job):
    if "files" in job:
        files = [BASE / relative for relative in job["files"]]
    else:
        files = []
        for folder in job["folders"]:
            folder = folder.strip()
            if not folder:
                continue
            steps = sorted(
                path
                for path in (BASE / folder).glob("*.py")
                if re.match(r"^_?\d+_", path.name)
            )
            if not steps:
                raise FileNotFoundError(
                    f"No numbered scripts found: {BASE / folder}"
                )
            files.extend(steps)

    if not files:
        raise RuntimeError("Analysis has no scripts")

    for path in files:
        if not path.is_file():
            raise FileNotFoundError(path)
        if not path.resolve().is_relative_to(BASE.resolve()):
            raise RuntimeError(f"Script outside mtec_high: {path}")

    return files


def dataset_ids(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if "dataset_id" not in (reader.fieldnames or []):
            raise ValueError(f"Missing dataset_id column: {path}")
        return {
            row["dataset_id"]
            for row in reader
            if row.get("dataset_id")
        }


def check_tsa_directional():
    """Reject missing or obviously stale legacy directional inputs."""
    folder = ROOT / "data/mtec_high/derived_data/gsea_tsa_excluded/gsea"

    for namespace in (
        "biological_process",
        "molecular_function",
        "cellular_component",
    ):
        stem = f"go_{namespace}_gsea"
        source = folder / "standard" / f"{stem}.csv"
        if not source.is_file():
            raise FileNotFoundError(source)

        expected_ids = dataset_ids(source)
        if not expected_ids:
            raise ValueError(f"No datasets in {source}")

        for direction in ("pos", "neg"):
            target = folder / direction / f"{stem}_{direction}.csv"
            if not target.is_file():
                raise FileNotFoundError(
                    f"Missing {target}. This branch currently has no "
                    "one-sided GSEA producer; fix that dependency first."
                )

            missing = expected_ids - dataset_ids(target)
            if missing:
                raise ValueError(
                    f"{target}: missing datasets {sorted(missing)}"
                )

            if target.stat().st_mtime_ns < source.stat().st_mtime_ns:
                raise RuntimeError(
                    f"Directional input predates refreshed GO results: "
                    f"{target}. Regenerate the TSA-excluded directional "
                    "inputs before retrying."
                )


def save_state(state):
    WORK.mkdir(parents=True, exist_ok=True)
    temporary = STATE_FILE.with_suffix(".tmp")
    temporary.write_text(
        json.dumps(state, indent=2),
        encoding="utf-8",
    )
    temporary.replace(STATE_FILE)


def descendants(names):
    affected = set(names)
    while True:
        expanded = affected | {
            name
            for name, job in JOBS.items()
            if affected.intersection(job["dependencies"])
        }
        if expanded == affected:
            return affected
        affected = expanded


def main():
    parser = argparse.ArgumentParser(
        description="Run mTEC-high analyses independently."
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--only", nargs="+", choices=list(JOBS))
    mode.add_argument("--retry", action="store_true")
    parser.add_argument("--list", action="store_true")
    args = parser.parse_args()

    if args.list:
        for name, job in JOBS.items():
            dependencies = ", ".join(job["dependencies"]) or "none"
            print(f"{name:26} dependencies: {dependencies}")
        return 0

    if not BASE.is_dir():
        parser.error("Save this script in the real analysis project root.")

    if args.retry or args.only:
        if STATE_FILE.is_file():
            state = json.loads(STATE_FILE.read_text(encoding="utf-8"))
        else:
            state = {}
        state = {
            name: state.get(name, {"status": "pending"})
            for name in JOBS
        }
    else:
        # A normal invocation starts a fresh full analysis run.
        state = {name: {"status": "pending"} for name in JOBS}

    if args.retry:
        selected = {
            name
            for name in JOBS
            if state[name]["status"] != "ok"
        }
    elif args.only:
        selected = set(args.only)
    else:
        selected = set(JOBS)

    if not selected:
        print("Nothing to retry.")
        return 0

    # Include prerequisites that have not previously completed.
    while True:
        additional = {
            dependency
            for name in selected
            for dependency in JOBS[name]["dependencies"]
            if state[dependency]["status"] != "ok"
        } - selected
        if not additional:
            break
        selected.update(additional)

    # Previously completed downstream outputs become stale when their
    # prerequisites are rerun. --retry will subsequently refresh them.
    for name in descendants(selected):
        state[name] = {"status": "pending"}
    save_state(state)

    run_id = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    logs = WORK / "logs" / run_id
    logs.mkdir(parents=True)

    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT) + (
        os.pathsep + env["PYTHONPATH"]
        if env.get("PYTHONPATH")
        else ""
    )
    env["PYTHONIOENCODING"] = "utf-8"
    env["PYTHONDONTWRITEBYTECODE"] = "1"

    # Some existing scripts write directly without creating their folders.
    # Create only known output directories; do not recreate obsolete paths.
    output_directories = (
        "genomics/window_sliding/source",
        "genomics/window_sliding/standardized",
        "genomics/curated_loci/source",
        "tsa_expression/aps1_known_antigens/standardized",
        "gsea_tsa_excluded/rodent_datasets_tsas_excluded",
        "gsea_tsa_excluded/gsea/standard/selected",
        "gsea_tsa_excluded/gsea/standard/nes_contrast",
    )
    for relative in output_directories:
        (ROOT / "data/mtec_high/derived_data" / relative).mkdir(
            parents=True, exist_ok=True
        )

    for name, job in JOBS.items():
        if name not in selected:
            continue

        unavailable = [
            dependency
            for dependency in job["dependencies"]
            if state[dependency]["status"] != "ok"
        ]
        if unavailable:
            message = f"Unsuccessful prerequisites: {unavailable}"
            state[name] = {"status": "blocked", "error": message}
            save_state(state)
            print(f"\nBLOCKED {name}: {message}", flush=True)
            continue

        log_file = logs / f"{name}.log"
        state[name] = {"status": "running", "log": str(log_file)}
        save_state(state)
        print(f"\nRUNNING {name}\nLog: {log_file}", flush=True)

        try:
            with log_file.open("w", encoding="utf-8") as log:
                try:
                    if job["guard"] == "tsa_directional":
                        check_tsa_directional()

                    for script in scripts_for(job):
                        module = ".".join(
                            script.relative_to(ROOT).with_suffix("").parts
                        )
                        print(f"  {script.relative_to(BASE)}", flush=True)
                        log.write(f"\n=== {module} ===\n")
                        log.flush()

                        subprocess.run(
                            [sys.executable, "-u", "-m", module],
                            cwd=ROOT,
                            env=env,
                            stdout=log,
                            stderr=subprocess.STDOUT,
                            check=True,
                        )
                except Exception as error:
                    log.write(f"\nMASTER RUNNER ERROR: {error}\n")
                    raise

        except KeyboardInterrupt:
            state[name]["status"] = "interrupted"
            save_state(state)
            print("\nInterrupted. Use --retry to resume.")
            return 130

        except Exception as error:
            state[name].update(status="failed", error=str(error))
            print(f"FAILED {name}: {error}", flush=True)

        else:
            state[name]["status"] = "ok"
            print(f"OK {name}", flush=True)

        save_state(state)

    print("\nAnalysis status:")
    for name in JOBS:
        print(f"  {state[name]['status']:12} {name}")

    incomplete = [
        name for name in JOBS if state[name]["status"] != "ok"
    ]
    print(f"\nLogs: {logs}")
    if incomplete:
        print(
            "Some analyses remain incomplete or stale. "
            "Fix reported problems, then run with --retry."
        )
    return 1 if incomplete else 0


if __name__ == "__main__":
    raise SystemExit(main())