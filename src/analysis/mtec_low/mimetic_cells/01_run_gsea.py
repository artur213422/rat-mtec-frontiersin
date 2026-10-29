from pathlib import Path

from src.analysis.mtec_high.isg_pathways_gsea.shared.gsea_helpers import attach_gene_symbols, prepare_ranked_genes
from src.common.dataset_catalog.api.get_datasets import get_mtec_low
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.excel import create_excel_table
from src.common.feratures.gsea.gsea import run_gsea
from src.common.feratures.gsea.io import gmt_to_gene_sets_dict
from src.common.feratures.gsea.models import AnalysedGseaResult, GseaInput, GseaResult
from src.common.feratures.gsea.save_result import save_gsea_result
from src.common.file_io import save_csv_to_excel
from src.common.paths import PROJECT_ROOT, find_file_in_project

##########################################################################################
# PREPARE DATA
##########################################################################################

output_folder: Path = PROJECT_ROOT / "data/mtec_low/derived_data/mimetic_cells"
output_file: Path = output_folder / "gsea_analysis.csv"

gmt_file: Path = find_file_in_project("mimetic_cells_ortholog.gmt")
gene_sets: dict[str, list[str]] = gmt_to_gene_sets_dict(gmt_file)

if not gene_sets:
    raise ValueError(f"No mimetic-cell gene sets found in {gmt_file}")


datasets: list[DatasetMetadata] = get_mtec_low()

##########################################################################################
# RUN GSEA
##########################################################################################

results: list[AnalysedGseaResult] = []
gene_symbols_by_dataset: dict[str, dict[str, str]] = {}

for dataset in datasets:
    print("Analysing dataset:", dataset.dataset_id)
    ranked_genes, id_to_symbol = prepare_ranked_genes(dataset)
    gene_symbols_by_dataset[dataset.dataset_id] = id_to_symbol

    gsea_result: list[GseaResult] = run_gsea(
        GseaInput(
            rnk=ranked_genes,
            gene_sets_dict=gene_sets,
            min_size=7,
            max_size=max(7, max(len(genes) for genes in gene_sets.values())),
        )
    )
    results.append(
        AnalysedGseaResult(
            gsea_result=gsea_result,
            dataset_id=dataset.dataset_id,
        )
    )



##########################################################################################
# SAVE
##########################################################################################

save_gsea_result(output_file, results)
attach_gene_symbols(output_file, gene_symbols_by_dataset)

excel_file: Path = save_csv_to_excel(output_file)
create_excel_table(excel_file)
print("Done:", output_file)
