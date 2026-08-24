import polars as pl
from pathlib import Path

from src.common.excel import create_excel_table
from src.common.feratures.gsea.gsea import run_gsea
from src.common.feratures.gsea.models import GseaResult, GseaInput, AnalysedGseaResult
from src.common.feratures.gsea.rnk import calculate_rnk_matrix, calculate_rnk_matrix_rat_special
from src.common.feratures.gsea.save_result import save_gsea_result
from src.common.file_io import save_csv_to_excel
from src.common.models.dataset_metadata.model import DatasetCatalog
from src.common.models.dataset_metadata.ops import build_dataset_catalog
from src.common.paths import resolve_dir_from_root, find_file_in_project
from src.common.types.enums import Species
from src.exploration.antigen_presentation_comparison.internal.methods import read_pathways_to_dict

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


pathways_file: Path = find_file_in_project("mapped_ifn_vs_other_ifn.csv")

datasets: DatasetCatalog = build_dataset_catalog()
output: Path = resolve_dir_from_root("data/analysis/isg_pathways_analysis/gse161664_human_lung_epithelia_cells/ifn_vs_other_ifn/result") / "gsea_analysis.csv"

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

result: list[AnalysedGseaResult] = []

for dataset in datasets.datasets:
    print("Analysing dataset:", dataset.dataset_id)
    tec_df: pl.DataFrame = dataset.get_standardized_df()

    if dataset.species == Species.RAT:
        rnk_df: pl.DataFrame = calculate_rnk_matrix_rat_special(tec_df)
    else:
        rnk_df: pl.DataFrame = calculate_rnk_matrix(tec_df)

    isg_pathways_df: pl.DataFrame = pl.read_csv(pathways_file)

    gsea_input: GseaInput = GseaInput(
        rnk=rnk_df,
        min_size=7,
        gene_sets_dict=read_pathways_to_dict(isg_pathways_df)
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
save_gsea_result(output, result)
print("Done")


excel_file: Path = save_csv_to_excel(output)
create_excel_table(excel_file)