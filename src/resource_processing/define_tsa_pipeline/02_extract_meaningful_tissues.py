from pathlib import Path

import polars as pl

from src.common.paths import try_find_file_in_project, resolve_dir_from_root
from src.resource_processing.define_tsa_pipeline.internal.tissue_extraction import extract_tissue_matrix, \
    extract_list_from_json

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

expression_matrix_file: Path = try_find_file_in_project("GSE10246_series_matrix.csv")

mouse_atlas_folder: Path = resolve_dir_from_root("resources/mouse_gene_atlas")

all_tissue_cols_file: Path = try_find_file_in_project("col_name_tissues.json")

picked_tissues_file: Path = mouse_atlas_folder / "standardized/picked_tissues.json"

if not picked_tissues_file.exists():
    raise FileNotFoundError(
        f"File {picked_tissues_file} does not exist. "
        "This file must be made manually and contain tissue names picked for TSA extraction."
    )

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

picked_tissues: list[str] = extract_list_from_json(picked_tissues_file)
all_tissue_cols: list[str] = extract_list_from_json(all_tissue_cols_file)

picked_matrix: pl.DataFrame = extract_tissue_matrix(
    expression_matrix_file=expression_matrix_file,
    picked_tissues=picked_tissues,
    all_tissue_cols=all_tissue_cols,
)

print(f"Input matrix: {expression_matrix_file}")
print(f"Picked tissues: {len(picked_tissues)}")
print(f"Output shape: {picked_matrix.shape}")

##########################################################################################
#                                           SAVE
##########################################################################################

output_file: Path = mouse_atlas_folder / "raw/GSE10246_picked_tissues_matrix.csv"

picked_matrix.write_csv(output_file)

print(f"Saved: {output_file}")