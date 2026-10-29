import json
from pathlib import Path

import polars as pl

from src.common.paths import find_file_in_project, resolve_dir_from_root
from src.resource_processing.define_tsa_pipeline.internal.tissue_extraction import extract_tissue_matrix, \
    extract_list_from_json

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

expression_matrix_file: Path = find_file_in_project("gse10246_with_replicate_averaged.csv")
expression_matrix: pl.DataFrame = pl.read_csv(expression_matrix_file)
out_file: Path = resolve_dir_from_root("resources/mouse_gene_atlas/intermediate") / "averaged_gse10246_after_excluding_samples.csv"


excluded_sample_group_file: Path = find_file_in_project("excluded_sample_groups.json")

with excluded_sample_group_file.open("r", encoding="utf-8") as file:
    grouped_tissue_map: list[str] = json.load(file)

##########################################################################################
# PROCESS DATA
##########################################################################################

excluded_sample_groups: list[str] = grouped_tissue_map
excluded_group_set: set[str] = set(excluded_sample_groups)

if "ID_REF" not in expression_matrix.columns:
    raise ValueError("Required ID_REF column is missing")

columns_to_drop: list[str] = [
    column
    for column in expression_matrix.columns
    if column in excluded_group_set
]

matched_groups: set[str] = set(columns_to_drop)

unmatched_groups: list[str] = sorted(excluded_group_set - matched_groups)

if unmatched_groups:
    raise ValueError(f"Excluded groups did not match matrix columns: {unmatched_groups}")

if not columns_to_drop:
    raise ValueError("No sample-group columns were selected for exclusion")

print(f"Dropping {len(columns_to_drop)} averaged sample groups:")

for column in columns_to_drop:
    print(f"- {column}")

expression_matrix = expression_matrix.drop(
    columns_to_drop
)

print(f"Remaining physiological groups: {expression_matrix.width - 1}")
##########################################################################################
#                                           SAVE
##########################################################################################

expression_matrix.write_csv(out_file)

print(f"Saved: {out_file}")