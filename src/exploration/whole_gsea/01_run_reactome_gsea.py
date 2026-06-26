

import polars as pl
from pathlib import Path

from src.common.excel import create_excel_table
from src.common.feratures.gsea.gsea import run_gsea
from src.common.feratures.gsea.io import gmt_to_gene_sets_dict
from src.common.feratures.gsea.models import GseaResult, GseaInput, AnalysedGseaResult
from src.common.feratures.gsea.rnk import calculate_rnk_matrix
from src.common.feratures.gsea.save_result import save_gsea_result
from src.common.io import save_csv_to_excel
from src.common.models.dataset_metadata.model import DatasetCatalog
from src.common.models.dataset_metadata.ops import build_dataset_catalog
from src.common.paths import resolve_dir_from_root, find_file_in_project
from src.exploration.antigen_presentation_comparison.internal.methods import read_pathways_to_dict

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


reactome_file: Path = find_file_in_project("reactome_pathways.gmt")

data_source: Path = resolve_dir_from_root("data/standardized")

datasets: DatasetCatalog = build_dataset_catalog()

reactome_gsea_output_file: Path = resolve_dir_from_root("data/analysis/whole_gsea") / "reactome_gsea.csv"

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

result: list[AnalysedGseaResult] = []

for dataset in datasets.datasets:
    file: Path = dataset.get_standardized_data_file(data_source)

    tec_df: pl.DataFrame = pl.read_csv(file)

    rnk_df: pl.DataFrame = calculate_rnk_matrix(tec_df)

    gsea_input: GseaInput = GseaInput(
        rnk=rnk_df,
        gene_sets_dict=gmt_to_gene_sets_dict(reactome_file)
    )

    gsea_result:  list[GseaResult] = run_gsea(gsea_input)

    result.append(
        AnalysedGseaResult(
            gsea_result=gsea_result,
            dataset_id=dataset.dataset_id
        )
    )

##########################################################################################
#                                       SAVE
##########################################################################################
save_gsea_result(reactome_gsea_output_file, result)
save_csv_to_excel(reactome_gsea_output_file)
create_excel_table(reactome_gsea_output_file)
print("Done")