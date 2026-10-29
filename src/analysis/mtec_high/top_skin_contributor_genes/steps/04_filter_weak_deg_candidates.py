from pathlib import Path

import polars as pl

from src.analysis.mtec_high.genomics.window_sliding.run_pipeline import PROJECT_ROOT
from src.analysis.mtec_high.top_skin_contributor_genes.config import CONTRIBUTOR_LOG2_FOLD_CHANGE_THRESHOLD, CONTRIBUTOR_P_VAL_THRESHOLD
from src.common.dataset_catalog.api.get_datasets import read_file_with_ortholog_id
from src.common.paths import find_file_in_project

# ------------------------------------------------------------------------------------------------
# Prepare Data
# ------------------------------------------------------------------------------------------------

output_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/top_skin_contributor_genes/standardized/top_skin_contributors.csv"
input_file: Path = find_file_in_project("03_cleaned_skin_contributor_genes.csv")

candidates: pl.DataFrame = read_file_with_ortholog_id(input_file)

# ------------------------------------------------------------------------------------------------
# Filter
# ------------------------------------------------------------------------------------------------

total_genes: int = candidates.height

candidates = candidates.filter((pl.col("p_val") < CONTRIBUTOR_P_VAL_THRESHOLD) & (pl.col("log2_fold_change").abs() >= CONTRIBUTOR_LOG2_FOLD_CHANGE_THRESHOLD))

remaining_genes: int = candidates.height
filtered_genes: int = total_genes - remaining_genes

print(f"Total candidates: {total_genes}")
print(f"Filtered out: {filtered_genes}")
print(f"Remaining candidates: {remaining_genes}")

# ------------------------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------------------------

candidates.write_csv(output_file)
