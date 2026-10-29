import polars as pl
from pathlib import Path

from src.analysis.mtec_high.tsa_expression_pattern.define_tsa_from_samsom_et_al.tissue_gsea.helpers.helpers import \
    make_tissue_gene_sets, get_mouse_symbol_ranking
from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.excel import create_excel_table
from src.common.feratures.gsea.gsea import run_gsea
from src.common.feratures.gsea.models import GseaResult, GseaInput, AnalysedGseaResult, GseaMode
from src.common.feratures.gsea.save_result import save_gsea_result
from src.common.file_io import save_csv_to_excel
from src.common.paths import resolve_dir_from_root, find_file_in_project

##########################################################################################
#                                       PREPARE DATA
##########################################################################################
source_file: Path = find_file_in_project("tsa_lists_with_rodent_expression.csv")
source_df: pl.DataFrame = pl.read_csv(source_file)
tissue_gene_sets: dict[str, list[str]] = make_tissue_gene_sets(source_df)

output_dir_file_neg: Path = resolve_dir_from_root("data/mtec_high/derived_data/tsa_expression/samsom_et_al_tsa_analysis/gsea_per_tissue") / "gsea_per_tissue_neg.csv"
output_dir_file_pos: Path = resolve_dir_from_root("data/mtec_high/derived_data/tsa_expression/samsom_et_al_tsa_analysis/gsea_per_tissue") / "gsea_per_tissue_pos.csv"

datasets: list[DatasetMetadata] = get_mtec_high()

result_pos: list[AnalysedGseaResult] = []
result_neg: list[AnalysedGseaResult] = []


##########################################################################################
#                                       RUN GSEA
##########################################################################################

for dataset in datasets:
    tec_df: pl.DataFrame = dataset.get_standardized_df()
    rnk_df: pl.DataFrame = get_mouse_symbol_ranking(dataset, tec_df)

    gsea_input = GseaInput(
        rnk=rnk_df,
        gene_sets_dict=tissue_gene_sets,
        min_size=6,
        permutation_num=1000,
        gsea_mode=GseaMode.Neg
    )

    gsea_result: list[GseaResult] = run_gsea(gsea_input)

    result_neg.append(
        AnalysedGseaResult(
            gsea_result=gsea_result,
            dataset_id=dataset.dataset_id,
        )
    )

save_gsea_result(output_dir_file_neg, result_neg)
excel_file: Path = save_csv_to_excel(output_dir_file_neg)
create_excel_table(excel_file)
print("Done")


for dataset in datasets:
    tec_df: pl.DataFrame = dataset.get_standardized_df()
    rnk_df: pl.DataFrame = get_mouse_symbol_ranking(dataset, tec_df)

    gsea_input = GseaInput(
        rnk=rnk_df,
        gene_sets_dict=tissue_gene_sets,
        min_size=6,
        permutation_num=1000,
        gsea_mode=GseaMode.Pos
    )

    gsea_result: list[GseaResult] = run_gsea(gsea_input)

    result_pos.append(
        AnalysedGseaResult(
            gsea_result=gsea_result,
            dataset_id=dataset.dataset_id,
        )
    )

save_gsea_result(output_dir_file_pos, result_pos)
excel_file: Path = save_csv_to_excel(output_dir_file_pos)
create_excel_table(excel_file)
print("Done")

##########################################################################################
#                                       SAVE
##########################################################################################


