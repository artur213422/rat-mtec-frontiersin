from pathlib import Path

import polars as pl

from src.common.dataframe_ops import extract_matrix
from src.common.feratures.analysers.deg_extractor import extract_top_genes
from src.common.feratures.analysers.deg_models import TopGenesResult
from src.common.models.dataset_metadata.model import DatasetCatalog
from src.common.models.dataset_metadata.ops import build_dataset_catalog
from src.common.paths import resolve_dir_from_root
from src.common.types.enums import Genotypes
from src.exploration.general_statistics.utils.analyse_degs import analyse_degs
from src.exploration.general_statistics.utils.calc_expression_statistics import calc_expression_statistics
from src.exploration.general_statistics.utils.calc_replicate_statistics import calc_replicate_statistics
from src.exploration.general_statistics.utils.dataset_statistrcs_json_parsing import save_statistics

from src.exploration.general_statistics.utils.types import DatasetStatistics, DegStatistics, ExpressionStatistics, \
    ReplicateStatistics

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


datasets: DatasetCatalog = build_dataset_catalog()

output: Path = resolve_dir_from_root("data/analysis/general_statistics")
save_file: Path = resolve_dir_from_root("data/analysis/general_statistics") / "dataset_statistics.json"

##########################################################################################
#                                       FILTER PATHWAYS
##########################################################################################

result: list[DatasetStatistics] = []

for dataset in datasets.datasets:
    print("Analyzing dataset: ", dataset.dataset_id)

    tec_df: pl.DataFrame = dataset.get_standardized_df()

    deg_statistics: DegStatistics = analyse_degs(tec_df)

    top_15: TopGenesResult = extract_top_genes(15, tec_df)


    ko_expression_statistics: ExpressionStatistics = calc_expression_statistics(
        extract_matrix(tec_df, Genotypes.KO)
    )

    ctrl_expression_statistics: ExpressionStatistics = calc_expression_statistics(
        extract_matrix(tec_df, Genotypes.CTRL)
    )

    combines_expression_statistics: ExpressionStatistics = calc_expression_statistics(
        extract_matrix(tec_df, Genotypes.ALL)
    )

    replicate_statistics: ReplicateStatistics = calc_replicate_statistics(
        extract_matrix(tec_df, Genotypes.CTRL),
        extract_matrix(tec_df, Genotypes.KO)
    )

    ##########################################################################################
    #                                       SAVE
    ##########################################################################################

    statistics: DatasetStatistics = DatasetStatistics(
        dataset_id=dataset.dataset_id,
        degs=deg_statistics,
        total_unique_genes=tec_df.height,

        ctrl=ctrl_expression_statistics,
        ko=ko_expression_statistics,
        combined=combines_expression_statistics,

        replicate_statistics=replicate_statistics,
        top_15=top_15
    )

    result.append(statistics)

save_statistics(result, save_file)




