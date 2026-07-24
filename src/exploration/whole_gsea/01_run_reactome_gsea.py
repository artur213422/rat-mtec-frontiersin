

import polars as pl
from pathlib import Path

from src.common.dataframe_ops import combine_datasets_with_nulls
from src.common.excel import create_excel_table
from src.common.feratures.gsea.gsea import run_gsea
from src.common.feratures.gsea.io import gmt_to_gene_sets_dict
from src.common.feratures.gsea.models import GseaResult, GseaInput, AnalysedGseaResult
from src.common.feratures.gsea.rnk import calculate_rnk_matrix
from src.common.feratures.gsea.save_result import save_singe_gsea_result
from src.common.file_io import save_csv_to_excel
from src.common.models.dataset_metadata.model import DatasetCatalog, DatasetMetadata
from src.common.models.dataset_metadata.ops import build_dataset_catalog
from src.common.paths import resolve_dir_from_root, find_file_in_project



reactome_file: Path = find_file_in_project("reactome_pathways.gmt")
datasets: DatasetCatalog = build_dataset_catalog()


##########################################################################################
#                                       RAT
##########################################################################################
####### PREPARE DATA

reactome_rat_gsea_output_file: Path = resolve_dir_from_root("data/analysis/whole_gsea") / "rat_reactome_gsea.csv"

rat_dataset: DatasetMetadata = datasets.get_rat_datasets()

rat_dataset_df: pl.DataFrame = rat_dataset.get_standardized_df()

####### PROCESS DATA

rat_result: list[AnalysedGseaResult] = []

rnk_df: pl.DataFrame = calculate_rnk_matrix(rat_dataset_df)

gsea_input: GseaInput = GseaInput(
    rnk=rnk_df,
    gene_sets_dict=gmt_to_gene_sets_dict(reactome_file)
)

rat_gsea_result: list[GseaResult] = run_gsea(gsea_input)



####### SAVE
save_singe_gsea_result(reactome_rat_gsea_output_file, rat_gsea_result)
rat_excel_file: Path = save_csv_to_excel(reactome_rat_gsea_output_file)
create_excel_table(rat_excel_file)
print("Done")

##########################################################################################
#                                       MOUSE
##########################################################################################
####### PREPARE DATA
reactome_mouse_gsea_output_file: Path = resolve_dir_from_root("data/analysis/whole_gsea") / "mouse_reactome_gsea.csv"

mouse_datasets: list[DatasetMetadata] = datasets.get_mouse_datasets()

combined_mouse_df: pl.DataFrame = combine_datasets_with_nulls(mouse_datasets)

####### PROCESS DATA

mouse_result: list[AnalysedGseaResult] = []

rnk_df: pl.DataFrame = calculate_rnk_matrix(rat_dataset_df)

gsea_input: GseaInput = GseaInput(
    rnk=rnk_df,
    gene_sets_dict=gmt_to_gene_sets_dict(reactome_file)
)

mouse_gsea_result: list[GseaResult] = run_gsea(gsea_input)



####### SAVE
save_singe_gsea_result(reactome_mouse_gsea_output_file, mouse_gsea_result)
mouse_excel_file: Path = save_csv_to_excel(reactome_rat_gsea_output_file)
create_excel_table(mouse_excel_file)
print("Done")