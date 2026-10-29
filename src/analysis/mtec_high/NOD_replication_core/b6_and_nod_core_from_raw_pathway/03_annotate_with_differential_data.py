from pathlib import Path

import polars as pl

from src.common.dataset_catalog.api.get_datasets import get_mouse_tec_high, get_rat_tec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.paths import PROJECT_ROOT
from src.common.types.enums import Species

# ------------------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------------------
input_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/NOD_replication_core/b6_and_nod_core_from_raw_pathway/source/02_pathway_genes_with_ortholog_id.csv"

output_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/NOD_replication_core/b6_and_nod_core_from_raw_pathway/source/03_annotated_pathways.csv"

mice_datasets: list[DatasetMetadata] = get_mouse_tec_high()
rat_dataset: DatasetMetadata = get_rat_tec_high()

# ------------------------------------------------------------------------------------------
# Annotate differential-expression statistics
# ------------------------------------------------------------------------------------------

statistic_columns: list[str] = [
    "log2_fold_change",
    "t_statistic",
    "p_val",
    "p_adj",
    "average_log2_expression",
]
output_columns: list[str] = [
    "dataset_id",
    "leading_gene",
    "ortholog_id",
    "gene_symbol",
    *statistic_columns,
]


def annotate_differential_data(input_file: Path, datasets: list[DatasetMetadata]) -> pl.DataFrame:
    genes_df = pl.read_csv(input_file, schema_overrides={"leading_gene": pl.String, "ortholog_id": pl.String},
    ).select("leading_gene", "ortholog_id").with_columns(
        pl.col("leading_gene").str.strip_chars().str.to_lowercase().alias("_gene_key"),
        pl.col("ortholog_id").str.strip_chars().replace("", None),
    )

    annotated_datasets: list[pl.DataFrame] = []
    for dataset in datasets:
        statistics = dataset.get_standardized_df().select(
            "gene_symbol",
            pl.col("ortholog_id").str.strip_chars().replace("", None),
            pl.col("gene_symbol").str.strip_chars().str.to_lowercase().alias("_gene_key"),
            *[pl.col(column).cast(pl.Float64) for column in statistic_columns],
        )

        if dataset.species is Species.MOUSE:
            # A matching source row takes priority even if its statistics are null.
            symbol_reference = statistics.drop("ortholog_id")
            matched = genes_df.join(symbol_reference, on="_gene_key", how="inner", validate="m:1")
            unmatched = genes_df.join(statistics.select("_gene_key"), on="_gene_key", how="anti")
        else:
            # Mouse input symbols must never be used to look up rat statistics.
            matched = None
            unmatched = genes_df

        # Preserve all source matches when an ortholog group has multiple genes.
        fallback = unmatched.join(
            statistics.drop("_gene_key"), on="ortholog_id", how="left"
        ).with_columns(pl.lit(dataset.dataset_id).alias("dataset_id")).select(output_columns)

        if matched is not None:
            matched = matched.with_columns(pl.lit(dataset.dataset_id).alias("dataset_id")).select(output_columns)
            fallback = pl.concat([matched, fallback])

        annotated_datasets.append(fallback)

    return pl.concat(annotated_datasets).sort("dataset_id", "leading_gene", "gene_symbol")


annotated_genes = annotate_differential_data(input_file, [*mice_datasets, rat_dataset])

# ------------------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------------------

output_file.parent.mkdir(parents=True, exist_ok=True)
annotated_genes.write_csv(output_file)
print(f"Saved {annotated_genes.height} annotated gene–dataset rows: {output_file}")
