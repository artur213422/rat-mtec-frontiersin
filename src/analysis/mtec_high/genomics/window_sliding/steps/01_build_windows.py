from pathlib import Path

import polars as pl

from src.analysis.mtec_high.genomics.window_sliding.config import WINDOW_MIN_SIZE, WINDOW_MAX_SIZE
from src.analysis.mtec_high.genomics.window_sliding.helpers.drop_invalid_chromosomes import drop_invalid_chromosomes
from src.analysis.mtec_high.genomics.window_sliding.helpers.enumerate_genomic_windows import enumerate_genomic_windows
from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.paths import PROJECT_ROOT

# -----------------------------------------------------------------------------------------------------------
# Prepare Data
# -----------------------------------------------------------------------------------------------------------

datasets: list[DatasetMetadata] = get_mtec_high()

output_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/genomics/window_sliding/source/01_windows.parquet"

# -----------------------------------------------------------------------------------------------------------
# Generate Windows
# -----------------------------------------------------------------------------------------------------------

windows: list[pl.DataFrame] = []

for dataset in datasets:
    tec_df: pl.DataFrame = dataset.get_standardized_df()

    tec_df = drop_invalid_chromosomes(tec_df, dataset.dataset_id)

    tec_df = enumerate_genomic_windows(tec_df=tec_df, dataset_id=dataset.dataset_id, min_size=WINDOW_MIN_SIZE, max_size=WINDOW_MAX_SIZE)
    windows.append(tec_df)

# -----------------------------------------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------------------------------------

result: pl.DataFrame = pl.concat(windows)

result.write_parquet(output_file)