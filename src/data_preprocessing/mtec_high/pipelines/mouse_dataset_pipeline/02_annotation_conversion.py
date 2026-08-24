from pathlib import Path

from src.common.feratures.gsm_handling.gsm_conversions import AnnotationResult, annotate_expression_matrix
from src.common.models.dataset_metadata.model import DatasetMetadata, DatasetCatalog
from src.common.models.dataset_metadata.ops import build_dataset_catalog
from src.common.paths import find_file_in_project, resolve_dir_from_root
from src.data_preprocessing.mtec_high.internal.mouse_array_annotations_tools import fix_GSE33878

##########################################################################################
#                                       PREPARE DATA
##########################################################################################
datasets: DatasetCatalog = build_dataset_catalog()
mouse_datasets: list[DatasetMetadata] = datasets.get_mouse_datasets()

output: Path = resolve_dir_from_root("data/raw/intermediate/annotated_expression_matrix/mouse")
raw_input: Path = resolve_dir_from_root("data/raw/intermediate/expression_matrix/mouse")

##########################################################################################
#                                       PROCESS FILE
##########################################################################################
for dataset in mouse_datasets:
    ##  Preparation

    raw_file_path: Path = find_file_in_project(f"{dataset.dataset_id}.txt", raw_input)

    matrix: AnnotationResult = annotate_expression_matrix(Path(dataset.annotation_table_file), raw_file_path)


    if dataset.dataset_id == "mouse_gse33878":
        matrix.table = fix_GSE33878(matrix.table)

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