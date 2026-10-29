import re
from pathlib import Path
import sys

# Support direct execution from any working directory.
PROJECT_ROOT = Path(__file__).resolve().parents[4]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
import polars as pl

from src.common.dataset_catalog.api.get_datasets import get_mtec_low
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.feratures.gene_pattern_extraction.methods import find_genes_by_pattern


##########################################################################################
#                                       PREPARE DATA
##########################################################################################


output_dir: Path = PROJECT_ROOT / "data/mtec_low/derived_data/chemokine_expression"
standardized_folder: Path = PROJECT_ROOT / "data/mtec_low/standardized"

datasets: list[DatasetMetadata] = get_mtec_low()

chemokine_patterns: list[re.Pattern[str]] = [
    re.compile(r"^Ccl\d+[a-z]?$", re.IGNORECASE),
    re.compile(r"^Cxcl\d+[a-z]?$", re.IGNORECASE),
    re.compile(r"^Xcl\d+[a-z]?$", re.IGNORECASE),
    re.compile(r"^Cx3cl\d+[a-z]?$", re.IGNORECASE),
]

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

chemokine_dfs: list[pl.DataFrame] = []

for dataset in datasets:
    tec_df: pl.DataFrame = dataset.get_standardized_df(root=standardized_folder)

    tec_df = find_genes_by_pattern(
        data=tec_df,
        gene_patterns=chemokine_patterns,
        gene_col="gene_symbol",
    )

    tec_df = tec_df.with_columns(
        pl.lit(dataset.dataset_id).alias("dataset_id")
    )

    chemokine_dfs.append(tec_df)

chemokines_df: pl.DataFrame = pl.concat(
    chemokine_dfs,
    how="diagonal",
)


##########################################################################################
#                                       SAVE
##########################################################################################

output_dir.mkdir(parents=True, exist_ok=True)
chemokines_df.write_csv(output_dir / "chemokines.csv")

