import json
from pathlib import Path
from typing import Any

from src.common.models.dataset_metadata.model import DatasetCatalog, DatasetMetadata
from src.common.models.dataset_metadata.ops import build_dataset_catalog
from src.common.paths import resolve_dir_from_root, find_file_in_project, collect_files_in_folder
from src.data_preprocessing.internal.R.run_R import run_rscript
from src.data_preprocessing.pipelines.mouse_dataset_pipeline.interials.helpers import gsm_sample_name_from_cel

##########################################################################################
#                                       PREPARE DATA
##########################################################################################
datasets: DatasetCatalog = build_dataset_catalog()
mouse_datasets: list[DatasetMetadata] = datasets.get_mouse_datasets()
raw_cel_folder: Path = resolve_dir_from_root("data/raw/cel_files/mouse")
output: Path = resolve_dir_from_root("data/raw/intermediate/expression_matrix/mouse")

r_exe = Path("C:/R/R/bin/Rscript.exe")
oligo_script: Path = find_file_in_project("cel_to_txt_oligo.R")
affy_script: Path = find_file_in_project("cel_to_txt_affy.R")
##########################################################################################
#                                       PROCESS DATA
##########################################################################################

for dataset in mouse_datasets:
    folder: Path = raw_cel_folder / Path(dataset.dataset_name)
    output_file: Path = output / f"{dataset.dataset_id}.txt"

    cel_files: list[Path] = collect_files_in_folder(folder)

    sample_names: list[str] = [
        gsm_sample_name_from_cel(file)
        for file in cel_files
    ]

    config: dict[str, Any] = {
        "cel_files": [str(path) for path in cel_files],
        "sample_names": sample_names,
        "output_file": str(output_file),
        "platform_package": dataset.platform_package
        ,"quote_output": True

    }

    if dataset.cel_protocol == "oligo":

        run_rscript(r_exe, oligo_script, config)
    elif dataset.cel_protocol == "affy":
         run_rscript(r_exe, affy_script, config)
    else:
        raise RuntimeError(f"Unknown protocol: {dataset.cel_protocol}")