from pathlib import Path
from typing import Any

from config.paths import RSCRIPT_EXE
from src.common.dataset_catalog.api.get_datasets import get_mouse_tec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.paths import find_file_in_project, resolve_dir_from_root
from src.data_preprocessing.mtec_high.internal.R.run_R import run_rscript


##########################################################################################
#                                       PREPARE DATA
##########################################################################################
mouse_datasets: list[DatasetMetadata] = get_mouse_tec_high()

source: Path = resolve_dir_from_root("data/mtec_high/source/intermediate/annotated_expression_matrix/mouse")

script: Path = find_file_in_project("attach_gene_location.R")
##########################################################################################
#                                       PROCESS FILE
##########################################################################################
for dataset in mouse_datasets:
    ##  Preparation

    source_file: Path = find_file_in_project(f"{dataset.dataset_id}.csv", source)

    config: dict[str, Any] = {
        "source_file": str(source_file),
        "output_file": str(source_file),
        "platform_package": dataset.platform_package,
    }

    run_rscript(RSCRIPT_EXE, script, config)

print("Done!")