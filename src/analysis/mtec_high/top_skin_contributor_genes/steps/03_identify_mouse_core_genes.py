from pathlib import Path

import polars as pl

from src.common.paths import PROJECT_ROOT
from src.analysis.mtec_high.top_skin_contributor_genes.helpers.identify_mouse_core import identify_mouse_core
from src.common.dataset_catalog.api.get_datasets import read_file_with_ortholog_id

# ------------------------------------------------------------------------------------------------
# Prepare Data
# ------------------------------------------------------------------------------------------------

output_folder: Path = PROJECT_ROOT / "data/mtec_high/derived_data/top_skin_contributor_genes/source"

input_file: Path = output_folder / "02_annotated_leading_genes.csv"
output_file: Path = output_folder / "03_cleaned_skin_contributor_genes.csv"

candidates: pl.DataFrame = read_file_with_ortholog_id(input_file)

# ------------------------------------------------------------------------------------------------
# Identify Mouse Core
# ------------------------------------------------------------------------------------------------

pathways: list[str] = candidates["pathway"].unique(maintain_order=True).to_list()
total_rows = candidates.height
rat_rows = candidates.filter(pl.col("dataset_id").str.contains("rat", literal=True)).height

if pathways:
    candidates = pl.concat([
        identify_mouse_core(candidates, pathway)
        for pathway in pathways
    ])

print(f"Input rows: {total_rows}")
print(f"Rat rows excluded: {rat_rows}")
print(f"Mouse rows below core threshold: {total_rows - rat_rows - candidates.height}")
print(f"Retained mouse core rows: {candidates.height}")

# ------------------------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------------------------

candidates.write_csv(output_file)


