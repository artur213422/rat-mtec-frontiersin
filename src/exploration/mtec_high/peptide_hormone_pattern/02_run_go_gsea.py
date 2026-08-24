import polars as pl
from pathlib import Path

from src.common.excel import create_excel_table
from src.common.feratures.gsea.gsea import run_gsea
from src.common.feratures.gsea.models import GseaResult, GseaInput, AnalysedGseaResult
from src.common.feratures.gsea.rnk import calculate_rnk_matrix
from src.common.feratures.gsea.save_result import save_gsea_result
from src.common.file_io import save_csv_to_excel
from src.common.models.dataset_metadata.model import DatasetCatalog
from src.common.models.dataset_metadata.ops import build_dataset_catalog
from src.common.paths import resolve_dir_from_root, find_file_in_project
from src.exploration.antigen_presentation_comparison.internal.methods import read_pathways_to_dict

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


go_pathways_file: Path = find_file_in_project("go_peptide_hormone_pathway_genes.csv")


datasets: DatasetCatalog = build_dataset_catalog()
output: Path = resolve_dir_from_root("data/analysis/peptide_hormone_pattern/results") / "go_gsea_analysis.csv"

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

result: list[AnalysedGseaResult] = []

for dataset in datasets.datasets:
    tec_df: pl.DataFrame = dataset.get_standardized_df()


    rnk_df: pl.DataFrame = calculate_rnk_matrix(tec_df)
    go_pathways_df: pl.DataFrame = pl.read_csv(go_pathways_file)

    gsea_input: GseaInput = GseaInput(
        rnk=rnk_df,
        gene_sets_dict=read_pathways_to_dict(go_pathways_df)
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