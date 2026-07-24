import re
from pathlib import Path

import polars as pl

from src.common.feratures.gene_pattern_extraction.methods import (
    find_genes_by_pattern,
)
from src.common.feratures.gene_pattern_extraction.save_csv import (
    write_report,
)
from src.common.models.dataset_metadata.model import DatasetCatalog
from src.common.models.dataset_metadata.ops import build_dataset_catalog
from src.common.paths import resolve_dir_from_root

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

source: Path = resolve_dir_from_root("data/standardized")

output_dir: Path = resolve_dir_from_root(
    "data/analysis/kallikrein_expression_pattern/kallikrein_lists"
)

report_file: Path = output_dir / "kallikrein_coverage_report.txt"

datasets: DatasetCatalog = build_dataset_catalog()

kallikrein_patterns: list[re.Pattern[str]] = [
    # Canonical tissue kallikrein / kallikrein-related peptidases:
    # Klk1, Klk2, ..., Klk15
    re.compile(r"^Klk(?:[1-9]|1[0-5])$", re.IGNORECASE),

    # Rodent kallikrein 1-related peptidase family:
    # Klk1b3, Klk1b4, Klk1b22, Klk1b27, etc.
    re.compile(r"^Klk1b\d+[a-z]?$", re.IGNORECASE),

    # Plasma kallikrein:
    # Klkb1
    re.compile(r"^Klkb1$", re.IGNORECASE),
]

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

kallikrein_dfs: list[pl.DataFrame] = []

for dataset in datasets.datasets:
    tec_df: pl.DataFrame = dataset.get_standardized_df()

    tec_df = find_genes_by_pattern(
        data=tec_df,
        gene_patterns=kallikrein_patterns,
        gene_col="gene_symbol",
    )

    tec_df = tec_df.with_columns(
        pl.lit(dataset.dataset_id).alias("dataset_id")
    )

    kallikrein_dfs.append(tec_df)

kallikreins_df: pl.DataFrame = pl.concat(
    kallikrein_dfs,
    how="diagonal",
)

##########################################################################################
#                                           SAVE
##########################################################################################

kallikreins_df.write_csv(
    output_dir / "kallikreins.csv"
)

write_report(
    kallikreins_df,
    report_file,
    "kallikrein",
)