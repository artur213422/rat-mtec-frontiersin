from pathlib import Path

import polars as pl

from src.common.paths import PROJECT_ROOT

# ------------------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------------------

input_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/NOD_replication_core/b6_and_nod_core_from_raw_pathway/source/04_annotated_pathway_cleaned.csv"
output_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/NOD_replication_core/b6_and_nod_core_from_raw_pathway/source/05_filtered_candidates.csv"

p_value_threshold: float = 0.05
fold_change_threshold: float = 0.5

# ------------------------------------------------------------------------------------------
# Keep orthologs with a DEG in at least one dataset and data in at least two
# ------------------------------------------------------------------------------------------

df = pl.read_csv(input_file, schema_overrides={"ortholog_id": pl.String})

is_deg = ((pl.col("p_val") <= p_value_threshold) & (pl.col("log2_fold_change").abs() >= fold_change_threshold)).fill_null(False)

gene_summary = df.group_by("ortholog_id").agg(is_deg.any().alias("has_deg"),  pl.col("dataset_id").drop_nulls().n_unique().alias("dataset_count"),
                                              )
retained_ids = gene_summary.filter(pl.col("has_deg") & (pl.col("dataset_count") > 1)).select("ortholog_id")

# Preserve all dataset rows for retained genes, including non-DEG rows.
filtered_df = df.join(retained_ids, on="ortholog_id", how="semi").sort("dataset_id", "ortholog_id", "leading_gene")

# ------------------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------------------

output_file.parent.mkdir(parents=True, exist_ok=True)
filtered_df.write_csv(output_file)

print(f"Removed {gene_summary.filter(~pl.col('has_deg')).height} orthologs with no DEGs")
print(f"Removed {gene_summary.filter(pl.col('has_deg') & (pl.col('dataset_count') <= 1)).height} remaining orthologs present in only one dataset")
print(f"Retained {retained_ids.height} orthologs and {filtered_df.height} rows")
print(f"Saved to {output_file}")
