import json
from pathlib import Path

import polars as pl

from src.common.paths import find_file_in_project, resolve_dir_from_root

##########################################################################################
#                                      PREPARE DATA
##########################################################################################
input_file: Path = find_file_in_project("GSE10246_affy_rma_log2_probe_by_sample.csv")
sample_titles_file: Path = find_file_in_project("sample_names.json")

output_file: Path = resolve_dir_from_root("resources/mouse_gene_atlas/intermediate") / "gse10246_with_sample_names.csv"

data: pl.DataFrame = pl.read_csv(input_file)

with sample_titles_file.open("r", encoding="utf-8") as file:
    grouped_tissue_map: dict[str, dict[str, str]] = json.load(file)

tissue_map: dict[str, str] = {
    gsm_accession: sample_name
    for group_samples in grouped_tissue_map.values()
    for gsm_accession, sample_name in group_samples.items()
}
##########################################################################################
#                                        ANNOTATE
##########################################################################################


sample_columns: list[str] = [
    column
    for column in data.columns
    if column != "ID_REF"
]
missing_mappings: list[str] = sorted(
    set(sample_columns) - set(tissue_map)
)

unexpected_mappings: list[str] = sorted(
    set(tissue_map) - set(sample_columns)
)

if missing_mappings:
    raise ValueError(f"Missing tissue-name mappings for columns: {missing_mappings}")

if unexpected_mappings:
    raise ValueError(f"JSON contains GSM accessions absent from the matrix: {unexpected_mappings}")

new_column_names: list[str] = [
    tissue_map[column]
    for column in sample_columns
]
if len(new_column_names) != len(set(new_column_names)):
    duplicated_names: list[str] = sorted(
        {
            name
            for name in new_column_names
            if new_column_names.count(name) > 1
        }
    )

    raise ValueError(f"Tissue mapping produces duplicate column names: {duplicated_names}")


data = data.rename(tissue_map)

##########################################################################################
#                                           SAVE
##########################################################################################
data.write_csv(output_file)