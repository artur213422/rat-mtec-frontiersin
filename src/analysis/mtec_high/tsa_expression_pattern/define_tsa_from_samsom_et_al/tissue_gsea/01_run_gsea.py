import polars as pl
from pathlib import Path

from src.analysis.mtec_high.tsa_expression_pattern.define_tsa_from_samsom_et_al.tissue_gsea.helpers.helpers import \
    get_mouse_symbol_ranking, make_tissue_gene_sets
from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.excel import create_excel_table
from src.common.feratures.gsea.gsea import run_gsea
from src.common.feratures.gsea.models import GseaResult, GseaInput, AnalysedGseaResult
from src.common.feratures.gsea.save_result import save_gsea_result
from src.common.file_io import save_csv_to_excel
from src.common.paths import resolve_dir_from_root, find_file_in_project


##########################################################################################
#                                       PREPARE DATA
##########################################################################################
source_file: Path = find_file_in_project("tsa_lists_with_rodent_expression.csv")
source_df: pl.DataFrame = pl.read_csv(source_file)
tissue_gene_sets: dict[str, list[str]] = make_tissue_gene_sets(source_df)

output_dir_file: Path = resolve_dir_from_root("data/mtec_high/derived_data/tsa_expression/samsom_et_al_tsa_analysis/gsea_per_tissue") / "gsea_per_tissue.csv"

datasets: list[DatasetMetadata] = get_mtec_high()

result: list[AnalysedGseaResult] = []


##########################################################################################
#                                       RUN GSEA
##########################################################################################
def cap_gene_sets_by_rank(gene_sets: dict[str, list[str]], rnk: pl.DataFrame, max_genes: int = 1000) -> dict[str, list[str]]:
    rank_by_gene = {
        gene.lower(): abs(float(score))
        for gene, score in rnk.select(
            "gene_symbol", "ranking_statistic"
        ).iter_rows()
    }

    capped = {}
    for tissue, genes in gene_sets.items():
        members = {
            gene.strip()
            for gene in genes
            if gene and gene.strip() and gene.lower() in rank_by_gene
        }
        capped[tissue] = sorted(
            members,
            key=lambda gene: (-rank_by_gene[gene.lower()], gene.lower()),
        )[:max_genes]

    return capped

for dataset in datasets:
    tec_df: pl.DataFrame = dataset.get_standardized_df()
    rnk_df: pl.DataFrame = get_mouse_symbol_ranking(dataset, tec_df)

    ranked_tissue_gene_sets = cap_gene_sets_by_rank(
        tissue_gene_sets, rnk_df, max_genes=500
    )

    gsea_input = GseaInput(
        rnk=rnk_df,
        gene_sets_dict=ranked_tissue_gene_sets,
        min_size=7,
        max_size=500,
        permutation_num=1000
    )

    gsea_result: list[GseaResult] = run_gsea(gsea_input)

    result.append(
        AnalysedGseaResult(
            gsea_result=gsea_result,
            dataset_id=dataset.dataset_id,
        )
    )



##########################################################################################
#                                       SAVE
##########################################################################################


save_gsea_result(output_dir_file, result)
excel_file: Path = save_csv_to_excel(output_dir_file)
create_excel_table(excel_file)
print("Done")
