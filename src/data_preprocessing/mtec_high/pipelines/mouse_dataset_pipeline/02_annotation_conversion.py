from pathlib import Path

from src.common.dataset_catalog.api.get_datasets import get_mouse_tec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.feratures.gsm_handling.gsm_conversions import AnnotationResult, annotate_expression_matrix
from src.common.paths import find_file_in_project, resolve_dir_from_root
from src.data_preprocessing.mtec_high.internal.mouse_array_annotations_tools import fix_weird_GPL6246_gene_symbols
from src.data_preprocessing.shared.try_rescue_chemokines import try_rescue_chemokines
from src.data_preprocessing.shared.try_rescue_kallikreins import try_rescue_kallikreins

##########################################################################################
#                                       PREPARE DATA
##########################################################################################
mouse_datasets: list[DatasetMetadata] = get_mouse_tec_high()

output: Path = resolve_dir_from_root("data/mtec_high/source/intermediate/annotated_expression_matrix/mouse")
raw_input: Path = resolve_dir_from_root("data/mtec_high/source/intermediate/expression_matrix/mouse")

##########################################################################################
#                                       PROCESS FILE
##########################################################################################
for dataset in mouse_datasets:
    ##  Preparation

    raw_file_path: Path = find_file_in_project(f"{dataset.dataset_id}.txt", raw_input)

    annotation_table_file: Path = find_file_in_project(dataset.annotation_table_file)

    matrix: AnnotationResult = annotate_expression_matrix(annotation_table_file, raw_file_path)

    if dataset.dataset_id == "mouse_gse33878_mtec_high":
        matrix.table = fix_weird_GPL6246_gene_symbols(matrix.table)

    matrix.table = try_rescue_chemokines(matrix.table)
    matrix.table = try_rescue_kallikreins(matrix.table)

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