from pathlib import Path

import polars as pl

from src.common.dataset_catalog.api.get_datasets import get_mouse_tec_low
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.paths import find_file_in_project, resolve_dir_from_root


##########################################################################################
#                                       PREPARE DATA
##########################################################################################

output_file: Path = resolve_dir_from_root("data/mtec_low/analysis/affy_vs_rna_seq") / "union_genes.csv"

output_file.parent.mkdir(parents=True, exist_ok=True,)

mouse_datasets: list[DatasetMetadata] = get_mouse_tec_low()

limma_gse2585: Path = find_file_in_project("GSE2585_limma.csv")

limma_gse65617: Path = find_file_in_project("GSE65617_limma.csv")

pydeseq2_gse205410: Path = find_file_in_project("mouse_gse205410_mtec_low.csv", resolve_dir_from_root("data/mtec_low/standardized"))



##########################################################################################
#                                       HELPERS
##########################################################################################
REQUIRED_COLUMNS: set[str] = {
    "gene_symbol",
    "log2_fold_change",
    "p_adj",
}


def read_effect_table(
    file: Path,
    *,
    dataset_id: str,
) -> pl.DataFrame:
    frame: pl.DataFrame = pl.read_csv(
        file,
        infer_schema_length=10_000,
        schema_overrides={
            "gene_symbol": pl.String,
            "log2_fold_change": pl.Float64,
            "p_adj": pl.Float64,
        },
        null_values=[
            "",
            "NA",
            "NaN",
            "nan",
            "NULL",
            "null",
        ],
    )

    missing_columns: set[str] = (
        REQUIRED_COLUMNS - set(frame.columns)
    )

    if missing_columns:
        raise ValueError(
            f"File '{file}' is missing columns: "
            f"{sorted(missing_columns)}"
        )

    result: pl.DataFrame = (
        frame
        .select(
            pl.col("gene_symbol")
            .cast(pl.String)
            .str.strip_chars(),

            pl.col("log2_fold_change")
            .cast(pl.Float64),

            pl.col("p_adj")
            .cast(pl.Float64),
        )
        .filter(
            pl.col("gene_symbol").is_not_null()
            & (pl.col("gene_symbol") != "")
        )
        .with_columns(
            pl.lit(dataset_id).alias("dataset_id")
        )
        .select(
            "gene_symbol",
            "log2_fold_change",
            "p_adj",
            "dataset_id",
        )
    )

    duplicated_genes: pl.DataFrame = (
        result
        .group_by("gene_symbol")
        .len()
        .filter(pl.col("len") > 1)
    )

    if duplicated_genes.height > 0:
        examples: list[str] = (
            duplicated_genes
            .get_column("gene_symbol")
            .head(10)
            .to_list()
        )

        raise ValueError(
            f"Dataset '{dataset_id}' contains duplicate genes. "
            f"Examples: {examples}"
        )

    return result



##########################################################################################
#                                       READ DATASETS
##########################################################################################

dataset_files: dict[str, Path] = {
    "mouse_gse2585_mtec_low": limma_gse2585,
    "mouse_gse65617_mtec_low": limma_gse65617,
    "mouse_gse205410_mtec_low": pydeseq2_gse205410
}

effect_tables: list[pl.DataFrame] = []

for dataset_id, file in dataset_files.items():
    print(f"Reading {dataset_id}: {file}")

    effect_tables.append(
        read_effect_table(
            file,
            dataset_id=dataset_id,
        )
    )


##########################################################################################
#                                 EXTRACT FOUR-WAY OVERLAP
##########################################################################################

all_effects_df: pl.DataFrame = pl.concat(
    effect_tables,
    how="vertical",
)

number_of_datasets: int = len(dataset_files)

valid_effects_df: pl.DataFrame = all_effects_df.filter(
    pl.col("log2_fold_change").is_not_null()
    & pl.col("log2_fold_change").is_not_nan()
    & pl.col("log2_fold_change").is_finite()
)

shared_genes_df: pl.DataFrame = (
    valid_effects_df
    .group_by("gene_symbol")
    .agg(
        pl.col("dataset_id")
        .n_unique()
        .alias("n_datasets")
    )
    .filter(
        pl.col("n_datasets") == number_of_datasets
    )
    .select("gene_symbol")
)

overlapped_df: pl.DataFrame = (
    valid_effects_df
    .join(
        shared_genes_df,
        on="gene_symbol",
        how="semi",
    )
    .sort(
        "gene_symbol",
        "dataset_id",
    )
)


##########################################################################################
#                                       VALIDATE
##########################################################################################

invalid_genes_df: pl.DataFrame = (
    overlapped_df
    .group_by("gene_symbol")
    .agg(
        pl.len().alias("n_rows"),
        pl.col("dataset_id")
        .n_unique()
        .alias("n_datasets"),

        pl.col("log2_fold_change")
        .is_null()
        .sum()
        .alias("n_null_effects"),

        pl.col("log2_fold_change")
        .is_nan()
        .sum()
        .alias("n_nan_effects"),

        pl.col("log2_fold_change")
        .is_infinite()
        .sum()
        .alias("n_infinite_effects"),
    )
    .filter(
        (pl.col("n_rows") != number_of_datasets)
        | (pl.col("n_datasets") != number_of_datasets)
        | (pl.col("n_null_effects") > 0)
        | (pl.col("n_nan_effects") > 0)
        | (pl.col("n_infinite_effects") > 0)
    )
)

if invalid_genes_df.height > 0:
    raise ValueError(
        "The overlap table contains incomplete genes:\n"
        f"{invalid_genes_df.head(10)}"
    )

expected_rows: int = (
    overlapped_df
    .get_column("gene_symbol")
    .n_unique()
    * number_of_datasets
)

if overlapped_df.height != expected_rows:
    raise ValueError(
        f"Expected {expected_rows:,} rows, "
        f"found {overlapped_df.height:,}."
    )


##########################################################################################
#                                       SAVE
##########################################################################################

overlapped_df.write_csv(output_file)

shared_gene_count: int = (
    overlapped_df
    .get_column("gene_symbol")
    .n_unique()
)

print(f"Saved: {output_file}")
print(f"Shared genes: {shared_gene_count:,}")
print(f"Rows: {overlapped_df.height:,}")
print(
    f"Expected rows: "
    f"{shared_gene_count * number_of_datasets:,}"
)

print(
    overlapped_df
    .group_by("dataset_id")
    .agg(
        pl.len().alias("n_rows"),
        pl.col("gene_symbol")
        .n_unique()
        .alias("n_genes"),

        pl.col("p_adj")
        .is_not_null()
        .sum()
        .alias("n_with_p_adj"),
    )
    .sort("dataset_id")
)