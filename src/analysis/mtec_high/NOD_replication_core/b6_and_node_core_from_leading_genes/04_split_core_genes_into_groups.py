from pathlib import Path

import polars as pl

from src.analysis.mtec_high.NOD_replication_core.b6_and_node_core_from_leading_genes.config import OTHERS_DATASETS, NOD_DATASETS
from src.common.dataset_catalog.api.get_datasets import read_file_with_ortholog_id
from src.common.paths import PROJECT_ROOT

# ------------------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------------------

standardized_folder: Path = PROJECT_ROOT / "data/mtec_high/derived_data/NOD_replication_core/b6_and_nod_core_from_leading_genes/standardized"
input_file: Path = standardized_folder / "nod_and_b6_core.csv"
output_file: Path = standardized_folder / "grouped_core_genes.csv"
summary_file: Path = standardized_folder / "gene_group_summary.csv"

# ------------------------------------------------------------------------------------------
# Calculate group medians; missing/non-finite effects do not contribute
# ------------------------------------------------------------------------------------------

df: pl.DataFrame = read_file_with_ortholog_id(input_file)

if df.select("leading_gene", "dataset_id").is_duplicated().any():
    raise ValueError("Expected one row per leading_gene and dataset_id.")

valid_effect = pl.col("log2_fold_change").is_not_null() & pl.col("log2_fold_change").is_finite()
nod_rows = pl.col("dataset_id").is_in(NOD_DATASETS) & valid_effect
other_rows = pl.col("dataset_id").is_in(OTHERS_DATASETS) & valid_effect

summary = df.group_by("leading_gene").agg(
    pl.col("log2_fold_change").filter(nod_rows).median().alias("nod_median_log2_fold_change"),
    pl.col("log2_fold_change").filter(other_rows).median().alias("others_median_log2_fold_change"),
    pl.col("dataset_id").filter(nod_rows).n_unique().alias("nod_dataset_count"),
    pl.col("dataset_id").filter(other_rows).n_unique().alias("others_dataset_count"),
)

nod = pl.col("nod_median_log2_fold_change")
others = pl.col("others_median_log2_fold_change")
summary = summary.with_columns(
    (others - nod).alias("contrast_score"),
    pl.when((nod > 0) & (others > 0)).then(pl.lit("both_up"))
    .when((nod < 0) & (others < 0)).then(pl.lit("both_down"))
    .when(((nod < 0) & (others > 0)) | ((nod > 0) & (others < 0)))
    .then(pl.lit("opposite_directions"))
    .otherwise(pl.lit("unclassified"))
    .alias("response_group"),
).sort("response_group", "contrast_score", "leading_gene", nulls_last=True)

# Preserve membership and all measurements, including passenger BALB/c rows.
grouped_df = df.join(summary, on="leading_gene", how="left", validate="m:1").sort(
    "response_group", "contrast_score", "leading_gene", "dataset_id", nulls_last=True
)

# ------------------------------------------------------------------------------------------
# Save annotated rows and a one-row-per-gene summary
# ------------------------------------------------------------------------------------------

grouped_df.write_csv(output_file)
summary.write_csv(summary_file)

for group, count in summary.group_by("response_group").len().sort("response_group").iter_rows():
    print(f"{group}: {count} genes")

print(f"Saved {grouped_df.height} annotated rows: {output_file}")
print(f"Saved {summary.height} gene summaries: {summary_file}")