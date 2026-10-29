import re
from pathlib import Path

# Support direct execution from any working directory.

import polars as pl

from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.feratures.gene_pattern_extraction.methods import find_genes_by_pattern
from src.common.paths import PROJECT_ROOT

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

output_dir: Path = PROJECT_ROOT / "data/mtec_high/derived_data/muc_expression"

standardized_folder: Path = PROJECT_ROOT / "data/mtec_high/standardized"

datasets: list[DatasetMetadata] = get_mtec_high()

muc_patterns: list[re.Pattern[str]] = [
    re.compile(r"^Muc\d+[a-z]?$", re.IGNORECASE),
    re.compile(r"^Mucl\d+$", re.IGNORECASE),
]

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

muc_dfs: list[pl.DataFrame] = []

for dataset in datasets:
    tec_df: pl.DataFrame = dataset.get_standardized_df(root=standardized_folder)

    tec_df = find_genes_by_pattern(
        data=tec_df,
        gene_patterns=muc_patterns,
        gene_col="gene_symbol",
    )

    tec_df = tec_df.with_columns(
        pl.lit(dataset.dataset_id).alias("dataset_id")
    )

    muc_dfs.append(tec_df)

muc_df: pl.DataFrame = pl.concat(
    muc_dfs,
    how="diagonal",
)

##########################################################################################
#                                           SAVE
##########################################################################################

muc_df.write_csv(output_dir / "muc.csv")
