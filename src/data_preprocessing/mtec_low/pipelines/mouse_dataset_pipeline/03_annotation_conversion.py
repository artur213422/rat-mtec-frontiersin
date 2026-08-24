from pathlib import Path

from src.common.dataset_catalog.api.get_datasets import get_mouse_tec_low
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.feratures.gsm_handling.gsm_conversions import AnnotationResult, annotate_expression_matrix
from src.common.paths import find_file_in_project, resolve_dir_from_root

##########################################################################################
#                                       PREPARE DATA
##########################################################################################
mouse_datasets: list[DatasetMetadata] = get_mouse_tec_low()

output: Path = resolve_dir_from_root("data/mtec_low/raw/intermediate/annotated_expression_matrix/mouse")
raw_input: Path = resolve_dir_from_root("data/mtec_low/raw/intermediate/expression_matrix/mouse")

##########################################################################################
#                                       PROCESS FILE
##########################################################################################
for dataset in mouse_datasets:
    ##  Preparation
    if dataset.dataset_name != "GSE2585":
        print("Skipping dataset:", dataset.dataset_id)
        continue

    raw_file_path: Path = find_file_in_project(f"{dataset.dataset_id}.txt", raw_input)
    annotation_table_file: Path = find_file_in_project(dataset.annotation_table_file)

    matrix: AnnotationResult = annotate_expression_matrix(annotation_table_file, raw_file_path)


    ##########################################################################################
    #                                           SAVE
    ##########################################################################################
    output_file: Path = output / f"{dataset.dataset_id}.csv"
    matrix.table.write_csv(output_file)

    conversion_result_file: Path = output / f"{dataset.dataset_id}_array_conversion_result.txt"

    conversion_result_file.write_text(
        f"Dataset: {dataset.dataset_name}\n"
        f"Missing gene count: {matrix.missing_gene_count}\n"
        f"Mapped gene count: {matrix.mapped_gene_count}\n",
        encoding="utf-8",
    )

print("Done!")