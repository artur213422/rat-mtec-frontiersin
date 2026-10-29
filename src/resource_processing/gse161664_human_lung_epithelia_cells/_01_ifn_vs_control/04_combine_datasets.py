from pathlib import Path
import polars as pl
from src.common.paths import resolve_dir_from_root

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


limma_result_folder: Path = resolve_dir_from_root("resources/interferon/gse161664_lung/intermediate/ifn_vs_control/limma_result")

output_file: Path = resolve_dir_from_root("resources/interferon/gse161664_lung/intermediate/ifn_vs_control") / "gse161664_lung_epithelia_ifn_vs_all.csv"

ifn_files = [
    limma_result_folder / f"{condition}_vs_ctrl.csv"
    for condition in ("ifna", "ifnb", "ifng", "ifnl")
]
##########################################################################################
#                                       PROCESS DATA
##########################################################################################

ifn_dfs: list[pl.DataFrame] = []

for file in ifn_files:
    ifn_df: pl.DataFrame = pl.read_csv(file)

    ifn_dfs.append(ifn_df)

combined_df: pl.DataFrame = pl.concat(ifn_dfs)
##########################################################################################
#                                           SAVE
##########################################################################################

combined_df.write_csv(output_file)