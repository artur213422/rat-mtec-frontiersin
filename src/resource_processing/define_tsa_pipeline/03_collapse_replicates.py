from pathlib import Path

import polars as pl

from src.common.paths import try_find_file_in_project, resolve_dir_from_root
from src.resource_processing.define_tsa_pipeline.internal.collapsing_tissue_cols import extract_tissue_replicates, \
    collapse_cols, log2_matrix
from src.resource_processing.define_tsa_pipeline.internal.tissue_extraction import extract_tissue_matrix, \
    extract_list_from_json

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

expression_matrix_file: Path = try_find_file_in_project("GSE10246_picked_tissues_matrix.csv")

output: Path = resolve_dir_from_root("resources/mouse_gene_atlas/standardized") / "standardized_tissues_matrix.csv"
picked_tissue_file = try_find_file_in_project("picked_tissues.json")
##########################################################################################
#                                       PROCESS DATA
##########################################################################################

tissue_list: list[str] = extract_list_from_json(picked_tissue_file)

matrix: pl.DataFrame = pl.read_csv(expression_matrix_file)


matrix = log2_matrix(matrix)


for tissue in tissue_list:
    tissue_pair: list[str] = extract_tissue_replicates(tissue, matrix)
    matrix = collapse_cols(matrix, tissue, tissue_pair)
##########################################################################################
#                                           SAVE
##########################################################################################

matrix.write_csv(output)



