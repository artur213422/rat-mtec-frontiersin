from pathlib import Path
import polars as pl
from src.common.paths import resolve_dir_from_root

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


source_folder: Path = resolve_dir_from_root("resources/interferon/gse28183_tec/intermediate/ifn_vs_all")


inf_file_folder: Path = source_folder / "limma_result"
output_file: Path = source_folder / "gse281833_tec_ifn_vs_all.csv"

ifn_files: list[Path] = [
    inf_file_folder / f"{condition}_vs_ctrl.csv"
    for condition in ("ifna", "ifng")
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
