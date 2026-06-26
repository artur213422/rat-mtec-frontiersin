from pathlib import Path

import polars as pl

from src.common.dataframe_ops import combine_datasets_with_nulls
from src.common.feratures.analysers.deg_extractor import extract_deg, extract_top_genes
from src.common.feratures.analysers.deg_models import TopGenesResult
from src.common.models.dataset_metadata.model import DatasetCatalog, DatasetMetadata
from src.common.models.dataset_metadata.ops import build_dataset_catalog
from src.common.paths import find_file_in_project, resolve_dir_from_root
from src.common.types.dataset import GeneReference

##########################################################################################
#                                       PREPARE DATA
##########################################################################################



data_source: Path = resolve_dir_from_root("data/standardized")

datasets: DatasetCatalog = build_dataset_catalog()

output_file: Path = resolve_dir_from_root("data/analysis/general_statistics") / "general_statistics.csv"

mouse_datasets: list[DatasetMetadata] = datasets.get_mouse_datasets()
rat_dataset: DatasetMetadata = datasets.get_rat_datasets()

rat_dataset_df = pl.read_csv(rat_dataset.get_standardized_data_file())
combined_mouse_df: pl.DataFrame = combine_datasets_with_nulls(mouse_datasets)


##########################################################################################
#                                       FILTER PATHWAYS
##########################################################################################
mouse_deg: list[GeneReference] = extract_deg(combined_mouse_df)
rat_deg: list[GeneReference] = extract_deg(rat_dataset_df)

rat_top_10: TopGenesResult = extract_top_genes(10, rat_dataset_df)
mouse_top_10: TopGenesResult = extract_top_genes(10, combined_mouse_df)

print("Mouse:", mouse_top_10.top_downregulated)
print("Rat:", rat_top_10.top_downregulated)
##########################################################################################
#                                       SAVE
##########################################################################################

