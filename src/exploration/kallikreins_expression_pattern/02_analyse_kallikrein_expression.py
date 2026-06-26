from pathlib import Path

import polars as pl

from src.common.analysers import genotype_comparison
from src.common.models.dataset_metadata.model import DatasetCatalog
from src.common.models.dataset_metadata.ops import build_dataset_catalog
from src.common.models.gene_expression_result.gene_expression_result import GeneExpressionResult
from src.common.paths import resolve_dir_from_root, resolve_file_from_root
from dataclasses import asdict

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

source: Path = resolve_dir_from_root("data/standardized")

all_kallikreins_file: Path = resolve_file_from_root("data/analysis/kallikrein_expression_pattern/kallikrein_lists/kallikreins_per_dataset.csv")

datasets: DatasetCatalog = build_dataset_catalog()

output_file: Path = resolve_dir_from_root("data/analysis/kallikrein_expression_pattern/results") / "gene_expression_results.csv"

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

results: list[GeneExpressionResult] = []

for dataset in datasets.datasets:
    dataset_file: Path = dataset.get_standardized_data_file(source)

    print(f"Reading {dataset_file.name}")
    dataset_id: str = dataset.dataset_id
    tec_df: pl.DataFrame = pl.read_csv(dataset_file)

    kallikrein_per_dataset: pl.DataFrame = pl.read_csv(all_kallikreins_file).filter(pl.col("dataset_id") == dataset_id)



    for row in kallikrein_per_dataset.iter_rows(named=True):
        ortholog_id: int = row["ortholog_id"]

        result: GeneExpressionResult = genotype_comparison(tec_df, ortholog_id, dataset_id)

        results.append(result)


    print("reading")



##########################################################################################
#                                       SAVE
##########################################################################################

results_df = pl.DataFrame(asdict(result) for result in results)

results_df.write_csv(output_file)