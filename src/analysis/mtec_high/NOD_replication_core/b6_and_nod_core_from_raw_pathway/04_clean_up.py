from pathlib import Path

import polars as pl

from src.common.paths import PROJECT_ROOT

# ------------------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------------------

input_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/NOD_replication_core/b6_and_nod_core_from_raw_pathway/source/03_annotated_pathways.csv"
output_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/NOD_replication_core/b6_and_nod_core_from_raw_pathway/source/04_annotated_pathway_cleaned.csv"


# ------------------------------------------------------------------------------------------
# Cleanup
# ------------------------------------------------------------------------------------------

df = pl.read_csv(
    input_file,
    schema_overrides={"ortholog_id": pl.String}
)

clean_df = df.filter(
    pl.col("ortholog_id").is_not_null()
    & (pl.col("ortholog_id").str.strip_chars() != "")
    & ~pl.col("ortholog_id").str.contains("///", literal=True)
    & pl.col("log2_fold_change").is_not_null()
)

# ------------------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------------------

clean_df.write_csv(output_file)

print(f"Removed {df.height - clean_df.height} rows")
print(f"Remaining {clean_df.height} rows")
print(f"Saved to {output_file}")
