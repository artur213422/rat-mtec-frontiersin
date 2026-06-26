from pathlib import Path

import polars as pl
from polars import DataFrame

from src.common.models.dataset_metadata.model import DatasetMetadata, DatasetCatalog
from src.common.models.dataset_metadata.ops import build_dataset_catalog
from src.common.paths import find_file_in_project, resolve_dir_from_root
from src.data_preprocessing.internal.mouse_array_annotations_tools import get_matrix, get_platform_id, \
    get_annotation_table, get_gene_symbol, fix_GSE33878

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

ens_conversion_table_file: Path = find_file_in_project("biomart_ens_rat_symbol_annotations.txt")
rgb_conversion_table_file: Path = find_file_in_project("rat_rgb_symbol_annotations.txt")

datasets: DatasetCatalog = build_dataset_catalog()
mouse_datasets: list[DatasetMetadata] = datasets.get_mouse_datasets()

output: Path = resolve_dir_from_root("data/intermediate/mouse")
raw_data_path: Path = resolve_dir_from_root("data/raw/mouse")

annotation_tables_folder: Path = resolve_dir_from_root("reference/array_platform_annotation_table")



##########################################################################################
#                                       PROCESS FILE
##########################################################################################
for dataset in mouse_datasets:
    ##  Preparation

    file_path: Path = find_file_in_project(dataset.file.name)
    matrix: DataFrame = get_matrix(file_path)

    platform_id: str = get_platform_id(file_path)
    print(f"Processing: {dataset.file.name} with platform id: {platform_id}")


    annot_table: DataFrame = get_annotation_table(annotation_tables_folder, platform_id)


    gene_symbols: list[str | None] = []

    missing_gene_count: int = 0
    mapped_gene_count: int = 0

    for probe_id in matrix["ID_REF"]:

        selected_cols: list[str] = ["gene_assignment", "Gene Symbol"]
        gene_name: str | None = get_gene_symbol(str(probe_id), annot_table, selected_cols) ## Or Gene Symbol

        if gene_name is None:
            missing_gene_count += 1
            gene_symbols.append(None)
            continue

        mapped_gene_count += 1
        gene_symbols.append(gene_name)

    matrix = matrix.with_columns(pl.Series("gene_symbol", gene_symbols))

    if platform_id == "GPL6246":
        matrix = fix_GSE33878(matrix)

    ##########################################################################################
    #                                           SAVE
    ##########################################################################################
    output_file: Path = output / f"{dataset.file.stem}.csv"
    matrix.write_csv(output_file)

    conversion_result_file: Path = output / f"{dataset.file.stem}_array_conversion_result.txt"

    conversion_result_file.write_text(
        f"Dataset: {dataset.dataset_name}\n"
        f"Missing gene count: {missing_gene_count}\n"
        f"Mapped gene count: {mapped_gene_count}\n",
        encoding="utf-8",
    )

print("Done!")