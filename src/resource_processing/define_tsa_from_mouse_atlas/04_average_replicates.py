import json
from pathlib import Path

import polars as pl

from src.common.paths import resolve_dir_from_root, find_file_in_project
from src.resource_processing.define_tsa_pipeline.internal.averaging_samples import average_sample

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

expression_matrix_file: Path = find_file_in_project("gse10246_with_sample_names.csv")

output: Path = resolve_dir_from_root("resources/mouse_gene_atlas/intermediate") / "gse10246_with_replicate_averaged.csv"

sample_titles_file: Path = find_file_in_project("sample_names.json")

with sample_titles_file.open("r", encoding="utf-8") as file:
    grouped_tissue_map: dict[str, dict[str, str]] = json.load(file)


expression_matrix: pl.DataFrame = pl.read_csv(expression_matrix_file)

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

for sample_group, group_samples in grouped_tissue_map.items():
    samples: list[str] = list(group_samples.values())

    expression_matrix = average_sample(matrix=expression_matrix, sample_columns=samples, output_column_name=sample_group)

##########################################################################################
#                                           SAVE
##########################################################################################

expression_matrix.write_csv(output)



