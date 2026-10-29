from pathlib import Path
import polars as pl

from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.paths import PROJECT_ROOT, find_file_in_project

##########################################################################################
#                                      PREPARE DATA
##########################################################################################

tsa_list: Path = find_file_in_project("mouse_tissue_atlas.csv")
tsa_list_df: pl.DataFrame = pl.read_csv(tsa_list)

tsa_ortholog_ids: list[str] = tsa_list_df.get_column("ortholog_id").cast(pl.String).unique().to_list()

output_dir: Path = PROJECT_ROOT / "data/mtec_high/derived_data/gsea_tsa_excluded/rodent_datasets_tsas_excluded"

datasets: list[DatasetMetadata] = get_mtec_high()

##########################################################################################
#                                        ANNOTATE
##########################################################################################

for dataset in datasets:
    output_file: Path = output_dir / f"{dataset.dataset_id}_tsa_excluded.csv"

    tec_df: pl.DataFrame = dataset.get_standardized_df()

    tec_df = tec_df.filter(
        ~pl.col("ortholog_id").is_in(tsa_ortholog_ids)
    )

    tec_df.write_csv(output_file)