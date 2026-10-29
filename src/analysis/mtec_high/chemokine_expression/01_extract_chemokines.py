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


output_dir: Path = PROJECT_ROOT / "data/mtec_high/derived_data/chemokine_expression"
standardized_folder: Path = PROJECT_ROOT / "data/mtec_high/standardized"

datasets: list[DatasetMetadata] = get_mtec_high()

chemokine_patterns: list[re.Pattern[str]] = [
    re.compile(r"^Ccl\d+[a-z]?(?:_family)?$", re.IGNORECASE),
    re.compile(r"^Cxcl\d+[a-z]?(?:_family)?$", re.IGNORECASE),
    re.compile(r"^Xcl\d+[a-z]?(?:_family)?$", re.IGNORECASE),
    re.compile(r"^Cx3cl\d+[a-z]?(?:_family)?$", re.IGNORECASE),
    re.compile(r"^Pf4$", re.IGNORECASE),
    re.compile(r"^Ppbp$", re.IGNORECASE),

]

def recover_chemokines(tec_df: pl.DataFrame) -> pl.DataFrame:
    symbol = pl.col("gene_symbol").str.to_lowercase()
    dataset = pl.col("dataset_id")

    return tec_df.with_columns(
        pl.when(
            ((symbol == "cxcl6") & dataset.str.starts_with("rat_"))
            | ((symbol == "cxcl5") & dataset.str.starts_with("mouse_"))
        )
        .then(pl.lit("Cxcl6_rat__Cxcl5_mouse"))
        .when(symbol == "pf4")
        .then(pl.lit("CXCL4"))
        .when(symbol == "ppbp")
        .then(pl.lit("CXCL7"))
        .otherwise(pl.col("gene_symbol"))
        .alias("gene_symbol")
    )

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
    tec_df = recover_chemokines(tec_df)
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

