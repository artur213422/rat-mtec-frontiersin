from pathlib import Path

from src.analysis.mtec_high.isg_pathways_gsea.shared.gsea_helpers import prepare_ranked_genes, prepare_isg_gene_sets, attach_gene_symbols
from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.excel import create_excel_table
from src.common.feratures.gsea.gsea import run_gsea
from src.common.feratures.gsea.models import AnalysedGseaResult, GseaInput, GseaResult
from src.common.feratures.gsea.save_result import save_gsea_result
from src.common.file_io import save_csv_to_excel
from src.common.paths import PROJECT_ROOT

##########################################################################################
# PREPARE DATA
##########################################################################################

input_output: dict[Path, Path] = {
    Path("resources/interferon/gse161664_lung/standardized/picked_isgs/gse161664_ifn_specific_induced_isgs.csv"):
        Path("data/mtec_high/derived_data/isg_pathways/gse161664_human_lung_epithelia_cells/gsea/ifn_specific_induced_isgs"),
    Path("resources/interferon/gse161664_lung/standardized/picked_isgs/gse161664_ifn_specific_intersection_isgs.csv"):
        Path("data/mtec_high/derived_data/isg_pathways/gse161664_human_lung_epithelia_cells/gsea/ifn_vs_ctrl_intersection"),
    Path("resources/interferon/gse161664_lung/standardized/picked_isgs/gse161664_ifn_vs_control_isgs.csv"):
        Path("data/mtec_high/derived_data/isg_pathways/gse161664_human_lung_epithelia_cells/gsea/ifn_vs_ctrl"),
}

datasets: list[DatasetMetadata] = get_mtec_high()

######################################################################################
# PROCESS DATA
######################################################################################

for relative_input_file, relative_out_dir in input_output.items():
    input_file: Path = PROJECT_ROOT / relative_input_file
    output: Path = PROJECT_ROOT / relative_out_dir / "gsea_analysis.csv"

    print("Preparing ISG gene sets:", input_file)
    isg_gene_sets: dict[str, list[str]] = prepare_isg_gene_sets(input_file)

    result: list[AnalysedGseaResult] = []
    gene_symbols_by_dataset: dict[str, dict[str, str]] = {}

    for dataset in datasets:
        print("Analysing dataset:", dataset.dataset_id)
        rnk_df, id_to_symbol = prepare_ranked_genes(dataset)
        gene_symbols_by_dataset[dataset.dataset_id] = id_to_symbol

        gsea_input = GseaInput(
            rnk=rnk_df,
            min_size=7,
            gene_sets_dict=isg_gene_sets,
            max_size=max(7, max(len(genes) for genes in isg_gene_sets.values())),
        )

        gsea_result: list[GseaResult] = run_gsea(gsea_input)
        result.append(
            AnalysedGseaResult(
                gsea_result=gsea_result,
                dataset_id=dataset.dataset_id,
            )
        )

    ##################################################################################
    # SAVE
    ##################################################################################

    save_gsea_result(output, result)
    attach_gene_symbols(output, gene_symbols_by_dataset)

    excel_file: Path = save_csv_to_excel(output)
    create_excel_table(excel_file)
    print("Done:", output)
