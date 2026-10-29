import polars as pl
from pathlib import Path
import sys

from src.common.paths import PROJECT_ROOT

# Support direct execution from any working directory.

from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.excel import create_excel_table
from src.common.feratures.gsea.gsea import run_gsea
from src.common.feratures.gsea.io import gmt_to_gene_sets_dict
from src.common.feratures.gsea.models import GseaResult, GseaInput, AnalysedGseaResult
from src.common.feratures.gsea.rnk import get_rnk_from_precalulated_source
from src.common.feratures.gsea.save_result import save_gsea_result
from src.common.file_io import save_csv_to_excel


##########################################################################################
#                                       PREPARE DATA
##########################################################################################

marker_gmt_file: Path = PROJECT_ROOT / "resources/hollander_2020_thymus_cell_markers/thymus_cell_markers.gmt"
output_dir_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/hollander_thymus_cell_types_gsea" / f"cell_type_gsea.csv"

datasets: list[DatasetMetadata] = get_mtec_high()

result: list[AnalysedGseaResult] = []


##########################################################################################
#                                       RUN GSEA
##########################################################################################


for dataset in datasets:

    tec_df: pl.DataFrame = dataset.get_standardized_df()

    rnk_df: pl.DataFrame = get_rnk_from_precalulated_source(tec_df)

    gsea_input: GseaInput = GseaInput(
        rnk=rnk_df,
        gene_sets_dict=gmt_to_gene_sets_dict(marker_gmt_file),
        max_size=600,

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
