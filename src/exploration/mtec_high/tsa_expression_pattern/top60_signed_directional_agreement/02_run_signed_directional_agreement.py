from pathlib import Path

import polars as pl

from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.paths import resolve_dir_from_root

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

source_folder: Path = resolve_dir_from_root("data/mtec_high/analysis/tsa_expression_pattern/signed_directional_agreement_against_rat")

source_file: Path = source_folder / "tsa_mtec_expression.csv"
out_file: Path = source_folder / "signed_directional_agreement.csv"

tsa_df: pl.DataFrame = pl.read_csv(source_file)

datasets: list[DatasetMetadata] = get_mtec_high()

##########################################################################################
#                                       SETTINGS
##########################################################################################

rat_dataset: str = "rat_mtec_high"
neutral_log2fc: float = 0.25

mouse_datasets: list[str] = [
    dataset.dataset_id
    for dataset in datasets
    if dataset.dataset_id != rat_dataset
]

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

rat_df: pl.DataFrame = (
    tsa_df
    .filter(
        (pl.col("dataset") == rat_dataset)
        & pl.col("found_in_dataset")
    )
    .select(
        "tissue_group",
        "gene_symbol",
        pl.col("log2_fold_change").alias("rat_log2_fold_change"),
    )
)

result_frames: list[pl.DataFrame] = []

for mouse_dataset in mouse_datasets:
    mouse_df: pl.DataFrame = (
        tsa_df
        .filter(
            (pl.col("dataset") == mouse_dataset)
            & pl.col("found_in_dataset")
        )
        .select(
            "tissue_group",
            "gene_symbol",
            pl.col("log2_fold_change").alias("mouse_log2_fold_change"),
        )
    )

    comparison_df: pl.DataFrame = (
        rat_df
        .join(
            mouse_df,
            on=["tissue_group", "gene_symbol"],
            how="inner",
        )
        .with_columns(
            pl.when(
                (
                    (pl.col("rat_log2_fold_change") > neutral_log2fc)
                    & (pl.col("mouse_log2_fold_change") > neutral_log2fc)
                )
                |
                (
                    (pl.col("rat_log2_fold_change") < -neutral_log2fc)
                    & (pl.col("mouse_log2_fold_change") < -neutral_log2fc)
                )
            )
            .then(1.0)
            .when(
                (
                    (pl.col("rat_log2_fold_change") > neutral_log2fc)
                    & (pl.col("mouse_log2_fold_change") < -neutral_log2fc)
                )
                |
                (
                    (pl.col("rat_log2_fold_change") < -neutral_log2fc)
                    & (pl.col("mouse_log2_fold_change") > neutral_log2fc)
                )
            )
            .then(-1.0)
            .otherwise(0.0)
            .alias("_directional_agreement")
        )
        .group_by("tissue_group")
        .agg(
            pl.col("_directional_agreement")
            .mean()
            .alias("agreement")
        )
        .with_columns(
            pl.lit(mouse_dataset).alias("dataset")
        )
        .select(
            "tissue_group",
            "dataset",
            "agreement",
        )
    )

    result_frames.append(comparison_df)

result_df: pl.DataFrame = pl.concat(result_frames).sort(["tissue_group", "dataset"])


##########################################################################################
#                                           SAVE
##########################################################################################

result_df.write_csv(out_file)