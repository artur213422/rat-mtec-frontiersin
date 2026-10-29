from pathlib import Path

import polars as pl

from src.analysis.mtec_high.genomics.window_sliding.columns import Windows
from src.analysis.mtec_high.genomics.window_sliding.helpers.score_genomic_windows import score_genomic_windows
from src.common.paths import PROJECT_ROOT, find_file_in_project

# -----------------------------------------------------------------------------------------------------------
# Prepare Data
# -----------------------------------------------------------------------------------------------------------

windows_file: Path = find_file_in_project("01_windows.parquet")
windows_df: pl.DataFrame = pl.read_parquet(windows_file)

if windows_df.schema != Windows.SCHEMA:
    raise ValueError(f"Invalid Windows schema.\n Expected: {Windows.SCHEMA}\n Got: {windows_df.schema}")

output_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/genomics/window_sliding/source/02_windows_score.parquet"

# -----------------------------------------------------------------------------------------------------------
# Score Windows
# -----------------------------------------------------------------------------------------------------------
datasets: list[str] = windows_df.get_column(Windows.DATASET_ID).unique().to_list()

scored_windows: list[pl.DataFrame] = []

for dataset in datasets:
    tec_df: pl.DataFrame = windows_df.filter(pl.col(Windows.DATASET_ID) == dataset)

    tec_df = score_genomic_windows(tec_df)
    scored_windows.append(tec_df)

# -----------------------------------------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------------------------------------

result: pl.DataFrame = pl.concat(scored_windows)

result.write_parquet(output_file)