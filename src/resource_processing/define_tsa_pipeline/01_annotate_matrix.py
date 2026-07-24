import json
from pathlib import Path

import polars as pl

from src.common.feratures.gsm_handling.gsm_conversions import AnnotationResult, annotate_expression_matrix, \
    get_platform_id, rename_expression_matrix_col, get_sample_names, clean_sample_names
from src.common.paths import try_find_file_in_project, resolve_dir_from_root

##########################################################################################
#                                       PREPARE DATA
##########################################################################################
raw_expression_matrix_file: Path = try_find_file_in_project("GSE10246_series_matrix.txt")

platform_id: str = "GPL1261"
annotation_table: Path = try_find_file_in_project("GPL1261-56135.txt")

output: Path = resolve_dir_from_root("resources/mouse_gene_atlas/raw")

actual_platform_id: str = get_platform_id(raw_expression_matrix_file)
if actual_platform_id != "GPL1261":
    raise ValueError(f"file platform id is {actual_platform_id}, but assumed platform_id")

output_json_file: Path = output / "col_name_tissues.json"
##########################################################################################
#                                       PROCESS DATA
##########################################################################################

annot_result: AnnotationResult = annotate_expression_matrix(annotation_table, raw_expression_matrix_file)

sample_names: list[str] = get_sample_names(raw_expression_matrix_file)

sample_names = clean_sample_names(sample_names)

matrix: pl.DataFrame = rename_expression_matrix_col(annot_result.table, sample_names)

matrix = matrix.drop("ID_REF")

##########################################################################################
#                                           SAVE
##########################################################################################

matrix.write_csv(output / "GSE10246_series_matrix.csv")

conversion_result_file: Path = output / f"annotation_result.txt"

conversion_result_file.write_text(
    f"Dataset: mouse atlas\n"
    f"Missing gene count: {annot_result.missing_gene_count}\n"
    f"Mapped gene count: {annot_result.mapped_gene_count}\n",
    encoding="utf-8",
)

with output_json_file.open("w", encoding="utf-8") as file:
    json.dump(sample_names, file, indent=4)
