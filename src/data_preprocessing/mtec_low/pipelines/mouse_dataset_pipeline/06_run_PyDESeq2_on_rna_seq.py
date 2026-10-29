from pathlib import Path

import polars as pl

from src.common.dataframe_ops import get_genotype_expression_values, get_genotype_col_names
from src.common.dataset_catalog.api.get_datasets import get_mouse_tec_low
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.feratures.PyDESeq2.run_pydeseq2 import run_pydeseq2, prepare_pydeseq2
from src.common.paths import resolve_dir_from_root, find_file_in_project
from src.common.types.enums import Genotypes

from src.data_preprocessing.mtec_high.internal.models.preprocessing_result.model import PreprocessingResult


##########################################################################################
#                                       PREPARE DATA
##########################################################################################

mouse_datasets: list[DatasetMetadata] = get_mouse_tec_low()


source_folder: Path = resolve_dir_from_root("data/mtec_low/standardized/mouse")
out_folder: Path = resolve_dir_from_root("data/mtec_low/standardized/mouse")


##########################################################################################
#                                       PROCESS FILE
##########################################################################################
for dataset in mouse_datasets:
    result: PreprocessingResult = PreprocessingResult(dataset.dataset_name)

    if dataset.dataset_name != "GSE205410":
        print(f"Skipping {dataset.dataset_name}")
        continue

    print(f"Running PyDESeq2 on {dataset.dataset_name}")

    tec_df_file: Path = find_file_in_project(f"{dataset.dataset_id}.csv", source_folder)
    tec_df: pl.DataFrame = pl.read_csv(tec_df_file)


    control_columns: list[str] = get_genotype_col_names(tec_df, Genotypes.CTRL)
    ko_columns: list[str] = get_genotype_col_names(tec_df, Genotypes.KO)

    tec_df = prepare_pydeseq2(counts=tec_df, control_columns=control_columns, ko_columns=ko_columns)
    tec_df = run_pydeseq2(counts=tec_df, control_columns=control_columns, ko_columns=ko_columns)


    ##########################################################################################
    #                                       SAVE
    ##########################################################################################

    tec_df.write_csv(tec_df_file)
