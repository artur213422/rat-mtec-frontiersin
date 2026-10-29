from pathlib import Path

import polars as pl

from src.common.dataset_catalog.api.get_datasets import get_mouse_tec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.feratures.limma.run_limma import LimmaDataType, run_limma, LimmaInput
from src.common.paths import resolve_dir_from_root, find_file_in_project


##########################################################################################
#                                       PREPARE DATA
##########################################################################################

mouse_datasets: list[DatasetMetadata] = get_mouse_tec_high()


source_folder: Path = resolve_dir_from_root("data/mtec_high/standardized/mouse")
out_folder: Path = resolve_dir_from_root("data/mtec_high/standardized/mouse")

r_script: Path = find_file_in_project("limma_script.R")

temporary_folder: Path = resolve_dir_from_root("data/temp/limma")

r_exe = Path("C:/R/R/bin/Rscript.exe")
##########################################################################################
#                                       PROCESS FILE
##########################################################################################


for dataset in mouse_datasets:
    tec_df: pl.DataFrame = dataset.get_standardized_df()

    output_file: Path = resolve_dir_from_root("data/mtec_high/standardized/mouse") / f"{dataset.dataset_id}.csv"

    limma_input: LimmaInput = LimmaInput(
        frame=tec_df,
        output_file=output_file,
        data_type="rma",
        merge=True
    )

    tec_df = run_limma(limma_input)

    ##########################################################################################
    #                                           SAVE
    ##########################################################################################


    tec_df.write_csv(output_file)
    print("Done")