
import polars as pl
from pathlib import Path

from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.excel import create_excel_table
from src.common.feratures.gsea.gsea import run_gsea
from src.common.feratures.gsea.models import GseaResult, GseaInput, AnalysedGseaResult
from src.common.feratures.gsea.rnk import calculate_rnk_matrix, calculate_rnk_matrix_rat_special, \
    get_rnk_from_precalulated_source
from src.common.feratures.gsea.save_result import save_gsea_result
from src.common.file_io import save_csv_to_excel
from src.common.paths import resolve_dir_from_root, find_file_in_project
from src.common.types.enums import Species
from src.exploration.mtec_high.antigen_presentation_comparison.internal.methods import read_pathways_to_dict

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


pathways_file: Path = find_file_in_project("gse281833_mapped_ifn_vs_all.csv")

datasets: list[DatasetMetadata] = get_mtec_high()
output: Path = resolve_dir_from_root("data/mtec_high/analysis/isg_pathways_analysis/gse281833_human_tec/ifn_vs_all_genes/result") / "gsea_analysis.csv"

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

result: list[AnalysedGseaResult] = []

for dataset in datasets:
    print("Analysing dataset:", dataset.dataset_id)
    tec_df: pl.DataFrame = dataset.get_standardized_df()

    rnk_df: pl.DataFrame = get_rnk_from_precalulated_source(tec_df)

    isg_pathways_df: pl.DataFrame = pl.read_csv(pathways_file)

    filtered_pathways_df: pl.DataFrame = pl.concat(
        [
            isg_pathways_df.filter(pl.col("pathway") != "ifng"),
            (
                isg_pathways_df
                .filter(pl.col("pathway") == "ifng")
                .sort("log2_fold_change", descending=True)
                .head(500)
            ),
        ],
        how="vertical",
    )

    gsea_input: GseaInput = GseaInput(
        rnk=rnk_df,
        min_size=7,
        gene_sets_dict=read_pathways_to_dict(filtered_pathways_df)
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