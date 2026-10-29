from pathlib import Path

import polars as pl

from src.common.excel import create_excel_table
from src.common.file_io import save_csv_to_excel
from src.common.paths import PROJECT_ROOT


input_dir = PROJECT_ROOT / "data/mtec_high/derived_data/gsea_tsa_excluded/gsea/standard"
output_dir = PROJECT_ROOT / "data/mtec_high/derived_data/gsea_tsa_excluded/gsea/standard/nes_contrast"

input_files = (
    "go_biological_process_gsea.csv",
    "go_molecular_function_gsea.csv",
    "go_cellular_component_gsea.csv",
)

for input_name in input_files:
    input_file = input_dir / input_name
    gsea_df = pl.read_csv(input_file)

    # Identify species from dataset_id.
    gsea_df = gsea_df.with_columns(
        pl.when(
            pl.col("dataset_id")
            .str.to_lowercase()
            .str.contains("rat")
        )
        .then(pl.lit("rat"))
        .when(
            pl.col("dataset_id")
            .str.to_lowercase()
            .str.contains("mouse|mice")
        )
        .then(pl.lit("mouse"))
        .otherwise(None)
        .alias("species")
    )

    # Count distinct datasets for each species.
    n_rat_datasets = (
        gsea_df
        .filter(pl.col("species") == "rat")
        .select("dataset_id")
        .n_unique()
    )

    n_mouse_datasets = (
        gsea_df
        .filter(pl.col("species") == "mouse")
        .select("dataset_id")
        .n_unique()
    )

    # Retain pathways where rat and ALL mice have opposite NES signs.
    retained_terms = (
        gsea_df
        .group_by("term")
        .agg(
            [
                pl.col("dataset_id")
                .filter(pl.col("species") == "rat")
                .n_unique()
                .alias("n_rat"),

                pl.col("dataset_id")
                .filter(pl.col("species") == "mouse")
                .n_unique()
                .alias("n_mouse"),

                (
                    pl.col("nes")
                    .filter(pl.col("species") == "rat")
                    > 0
                )
                .all()
                .alias("all_rat_positive"),

                (
                    pl.col("nes")
                    .filter(pl.col("species") == "rat")
                    < 0
                )
                .all()
                .alias("all_rat_negative"),

                (
                    pl.col("nes")
                    .filter(pl.col("species") == "mouse")
                    > 0
                )
                .all()
                .alias("all_mouse_positive"),

                (
                    pl.col("nes")
                    .filter(pl.col("species") == "mouse")
                    < 0
                )
                .all()
                .alias("all_mouse_negative"),
            ]
        )
        .filter(
            (pl.col("n_rat") == n_rat_datasets)
            & (pl.col("n_mouse") == n_mouse_datasets)
            & (
                (
                    pl.col("all_rat_positive")
                    & pl.col("all_mouse_negative")
                )
                |
                (
                    pl.col("all_rat_negative")
                    & pl.col("all_mouse_positive")
                )
            )
        )
        .select("term")
    )

    # Keep all dataset rows for retained pathways.
    filtered_df = (
        gsea_df
        .join(retained_terms, on="term", how="semi")
        .drop("species")
        .sort(["term", "dataset_id"])
    )

    output_file = output_dir / input_name.replace(
        ".csv",
        "_rat_vs_all_mice_nes_contrast.csv",
    )

    filtered_df.write_csv(output_file)

    output_excel: Path = save_csv_to_excel(output_file)
    create_excel_table(output_excel)

    print(f"{input_name}: retained {retained_terms.height} pathways; saved {filtered_df.height} rows to {output_file}")