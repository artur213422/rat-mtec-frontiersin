import polars as pl
from pathlib import Path
import sys

# Support direct execution from any working directory.
PROJECT_ROOT = Path(__file__).resolve().parents[4]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

ANALYSIS_FOLDER = PROJECT_ROOT / "data/mtec_high/derived_data/antigen_presentation_gsea"
STANDARDIZED_FOLDER = PROJECT_ROOT / "data/mtec_high/standardized"

from src.common.excel import create_excel_table
from src.common.feratures.gsea.gsea import run_gsea
from src.common.feratures.gsea.models import GseaResult, GseaInput, AnalysedGseaResult
from src.common.feratures.gsea.rnk import calculate_rnk_matrix, calculate_rnk_matrix_rat_special
from src.common.feratures.gsea.save_result import save_gsea_result
from src.common.file_io import save_csv_to_excel
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.types.enums import Species
from src.analysis.mtec_high.antigen_presentation_gsea.internal.methods import read_pathways_to_dict

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


pathways_file: Path = (ANALYSIS_FOLDER / "pathways/reactome_antigen_presentation_pathway_genes.csv")

datasets: list[DatasetMetadata] = get_mtec_high()
output: Path = (ANALYSIS_FOLDER / "gsea") / "reactome_gsea_analysis.csv"

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

result: list[AnalysedGseaResult] = []

for dataset in datasets:
    tec_df: pl.DataFrame = dataset.get_standardized_df(root=STANDARDIZED_FOLDER)

    if dataset.species == Species.RAT:
        rnk_df: pl.DataFrame = calculate_rnk_matrix_rat_special(tec_df)
    else:
        rnk_df: pl.DataFrame = calculate_rnk_matrix(tec_df)



    reactome_pathways_df: pl.DataFrame = pl.read_csv(pathways_file)

    gsea_input: GseaInput = GseaInput(
        rnk=rnk_df,
        min_size=7,
        gene_sets_dict=read_pathways_to_dict(reactome_pathways_df)
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