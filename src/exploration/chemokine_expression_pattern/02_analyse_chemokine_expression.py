from pathlib import Path

import polars as pl

from src.common.analysers import genotype_comparison
from src.common.paths import (
    resolve_dir_from_root,
    resolve_file_from_root,
)

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

kallikreins_file: Path = resolve_file_from_root(
    "data/analysis/kallikrein_expression_pattern/"
    "kallikrein_lists/kallikreins.csv"
)

kallikreins_df: pl.DataFrame = pl.read_csv(
    kallikreins_file,
    schema_overrides={
        "ortholog_id": pl.String,
    },
)

output_file: Path = (
    resolve_dir_from_root(
        "data/analysis/kallikrein_expression_pattern/"
        "kallikrein_lists"
    )
    / "analysed_kallikreins.csv"
)

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

kallikreins_df = genotype_comparison(kallikreins_df)

##########################################################################################
#                                           SAVE
##########################################################################################

kallikreins_df.write_csv(output_file)