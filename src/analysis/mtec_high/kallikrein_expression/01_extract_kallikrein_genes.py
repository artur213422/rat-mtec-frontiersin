import re
from pathlib import Path
import sys

# Support direct execution from any working directory.
PROJECT_ROOT = Path(__file__).resolve().parents[4]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
import polars as pl

from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.feratures.gene_pattern_extraction.methods import find_genes_by_pattern

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


output_dir: Path = (PROJECT_ROOT / "data/mtec_high/derived_data/kallikrein_expression")

standardized_folder: Path = PROJECT_ROOT / "data/mtec_high/standardized"

datasets: list[DatasetMetadata] = get_mtec_high()

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

for dataset in datasets:
    tec_df: pl.DataFrame = dataset.get_standardized_df(root=standardized_folder)

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

output_dir.mkdir(parents=True, exist_ok=True)
kallikreins_df.write_csv(output_dir / "kallikreins.csv")
