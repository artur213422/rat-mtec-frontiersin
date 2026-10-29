from pathlib import Path

import polars as pl

from src.common.paths import PROJECT_ROOT
from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata

# ------------------------------------------------------------------------------------------------
# Prepare Data
# ------------------------------------------------------------------------------------------------

dataset_metadata: list[DatasetMetadata] = get_mtec_high()

datasets: dict[str, pl.DataFrame] = {
    metadata.dataset_id: metadata.get_standardized_df()
    for metadata in dataset_metadata
}

output_folder: Path = PROJECT_ROOT / "data/mtec_high/derived_data/top_skin_contributor_genes/source"

input_file: Path = output_folder / "01_skin_leading_genes.csv"
output_file: Path = output_folder / "02_annotated_leading_genes.csv"

candidates: pl.DataFrame = pl.read_csv(input_file)


# ------------------------------------------------------------------------------------------------
# Annotate
# ------------------------------------------------------------------------------------------------
selected_columns: list[str] = [
    "gene_symbol",
    "chromosome",
    "strand",
    "start",
    "stop",
    "log2_fold_change",
    "t_statistic",
    "p_val",
    "p_adj",
    "average_log2_expression",
]

de_data: pl.DataFrame = pl.concat([
    dataset.select(selected_columns).with_columns(pl.lit(dataset_id).alias("dataset_id"))
    for dataset_id, dataset in datasets.items()
])

candidates = candidates.rename({"leading_gene": "gene_symbol"})

gene_key = pl.col("gene_symbol").str.strip_chars().str.to_lowercase().alias("_gene_key")
candidates = candidates.with_columns(gene_key)
de_data = de_data.with_columns(gene_key)

join_columns: list[str] = ["dataset_id", "_gene_key"]

# Check for unmatched genes
unmatched: pl.DataFrame = candidates.select(join_columns).unique().join(
    de_data.select(join_columns),
    on=join_columns,
    how="anti",
)

if not unmatched.is_empty():
    raise ValueError(f"Found {unmatched.height} genes without DE data:\n{unmatched}")

# Annotate candidates
candidates = candidates.drop("gene_symbol").join(
    de_data,
    on=join_columns,
    how="left",
    validate="m:1",
).drop("_gene_key")

# ------------------------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------------------------

candidates.write_csv(output_file)
