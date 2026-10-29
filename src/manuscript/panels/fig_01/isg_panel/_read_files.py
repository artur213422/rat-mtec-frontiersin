import csv
import math
from pathlib import Path

import polars as pl

from src.common.paths import PROJECT_ROOT
from src.manuscript.shared.config.dataset_labels import DATASET_LABELS
INPUT_FOLDER = PROJECT_ROOT / "data/mtec_high/derived_data/isg_pathways/gse281833_human_tec"

IFN_VS_CTRL_FILE: Path = INPUT_FOLDER / "ifn_vs_ctrl_gsea/gsea_analysis.csv"
INDUCED_ISG_FILE: Path = INPUT_FOLDER / "ifn_specific_induced_isgs_gsea/gsea_analysis.csv"


def _read_ifnb() -> pl.DataFrame:
    rows: list[dict[str, object]] = []

    with IFN_VS_CTRL_FILE.open(newline="", encoding="utf-8-sig",) as handle:

        reader = csv.DictReader(handle)

        for row in reader:

            if row["term"] != "ifng":
                continue

            if row["dataset_id"] not in DATASET_LABELS:
                continue

            rows.append(
                {
                    "dataset": row["dataset_id"],
                    "nes": float(row["nes"]),
                    "fdr": (
                        float(row["fdr_q_value"])
                        if row["fdr_q_value"].strip()
                        and row["fdr_q_value"].upper() != "NA"
                        else math.nan
                    ),
                }
            )

    return pl.DataFrame(rows)


def _read_ifna() -> pl.DataFrame:
    rows: list[dict[str, object]] = []

    with IFN_VS_CTRL_FILE.open(newline="", encoding="utf-8-sig") as handle:

        reader = csv.DictReader(handle)

        for row in reader:

            if row["term"] != "ifna":
                continue

            if row["dataset_id"] not in DATASET_LABELS:
                continue

            rows.append(
                {
                    "dataset": row["dataset_id"],
                    "nes": float(row["nes"]),
                    "fdr": (
                        float(row["fdr_q_value"])
                        if row["fdr_q_value"].strip()
                        and row["fdr_q_value"].upper() != "NA"
                        else math.nan
                    ),
                }
            )

    return pl.DataFrame(rows)