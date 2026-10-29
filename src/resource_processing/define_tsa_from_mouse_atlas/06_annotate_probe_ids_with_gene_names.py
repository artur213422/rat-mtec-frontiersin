import json
from pathlib import Path

import polars as pl

from src.common.feratures.gsm_handling.gsm_conversions import AnnotationResult, annotate_expression_matrix, \
    get_platform_id, rename_expression_matrix_col, get_sample_names, clean_sample_names
from src.common.paths import find_file_in_project, resolve_dir_from_root

##########################################################################################
#                                       PREPARE DATA
##########################################################################################
matrix_file: Path = find_file_in_project("averaged_gse10246_after_excluding_samples.csv")

platform_id: str = "GPL1261"
annotation_table: Path = find_file_in_project("GPL1261-56135.txt")

output: Path = resolve_dir_from_root("resources/mouse_gene_atlas/raw")


out_file: Path = resolve_dir_from_root("resources/mouse_gene_atlas/standardized") / "gse10246_atlas.csv"

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

annot_result: AnnotationResult = annotate_expression_matrix(annotation_table, matrix_file)

##########################################################################################
#                                           SAVE
##########################################################################################

annot_result.table.write_csv(out_file)

conversion_result_file: Path = output / f"annotation_result.txt"

conversion_result_file.write_text(
    f"Dataset: mouse atlas\n"
    f"Missing gene count: {annot_result.missing_gene_count}\n"
    f"Mapped gene count: {annot_result.mapped_gene_count}\n",
    encoding="utf-8",
)

