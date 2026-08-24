import polars as pl
from pathlib import Path

from src.common.dataset_catalog.api.get_datasets import  get_mtec_low
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.excel import create_excel_table
from src.common.feratures.gsea.gsea import run_gsea
from src.common.feratures.gsea.io import gmt_to_gene_sets_dict
from src.common.feratures.gsea.models import GseaResult, GseaInput, AnalysedGseaResult
from src.common.feratures.gsea.rnk import calculate_rnk_matrix, get_rnk_from_precalulated_source
from src.common.feratures.gsea.save_result import save_singe_gsea_result, save_gsea_result
from src.common.file_io import save_csv_to_excel
from src.common.paths import resolve_dir_from_root, find_file_in_project


##########################################################################################
#                                       PREPARE DATA
##########################################################################################

reactome_file: Path = find_file_in_project("thymus_cell_markers_hollander_2020.gmt")
output_dir_file: Path = resolve_dir_from_root("data/mtec_low/analysis/post_aire_analysis") / f"cell_type_gsea.csv"

datasets: list[DatasetMetadata] = get_mtec_low()

result: list[AnalysedGseaResult] = []


##########################################################################################
#                                       RUN GSEA
##########################################################################################


for dataset in datasets:

    tec_df: pl.DataFrame = dataset.get_standardized_df()

    rnk_df: pl.DataFrame = get_rnk_from_precalulated_source(tec_df)

    gsea_input: GseaInput = GseaInput(
        rnk=rnk_df,
        gene_sets_dict=gmt_to_gene_sets_dict(reactome_file)
    )

    gsea_result: list[GseaResult] = run_gsea(gsea_input)

    result.append(AnalysedGseaResult(
        gsea_result=gsea_result,
        dataset_id=dataset.dataset_id,
    ))


##########################################################################################
#                                       SAVE
##########################################################################################


save_gsea_result(output_dir_file, result)
excel_file: Path = save_csv_to_excel(output_dir_file)
create_excel_table(excel_file)
print("Done")
