from pathlib import Path

import polars as pl

from src.analysis.mtec_high.genomics.window_sliding.columns import ScoredWindows, Windows
from src.analysis.mtec_high.genomics.window_sliding.helpers.call_clusters import select_non_overlapping_clusters
from src.common.paths import PROJECT_ROOT, find_file_in_project

# -----------------------------------------------------------------------------------------------------------
# Prepare Data
# -----------------------------------------------------------------------------------------------------------

windows_file: Path = find_file_in_project("02_windows_score.parquet")
windows_df: pl.DataFrame = pl.read_parquet(windows_file)

if windows_df.schema != ScoredWindows.SCHEMA:
    raise ValueError(f"Invalid ScoredWindows schema.\n Expected: {ScoredWindows.SCHEMA}\n Got: {windows_df.schema}")

output_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/genomics/window_sliding/source/03_clusters.parquet"

# -----------------------------------------------------------------------------------------------------------
# Score Windows
# -----------------------------------------------------------------------------------------------------------
datasets: list[str] = windows_df.get_column(Windows.DATASET_ID).unique().to_list()

clustered: list[pl.DataFrame] = []

for dataset in datasets:
    tec_df: pl.DataFrame = windows_df.filter(pl.col(Windows.DATASET_ID) == dataset)

    tec_df = select_non_overlapping_clusters(tec_df)
    clustered.append(tec_df)

# -----------------------------------------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------------------------------------

result: pl.DataFrame = pl.concat(clustered)

result.write_parquet(output_file)