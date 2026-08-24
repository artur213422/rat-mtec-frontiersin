from pathlib import Path

import polars as pl

from src.common.dataset_catalog.api.get_datasets import (
    get_mouse_tec_low,
)
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.paths import resolve_dir_from_root


def standardize_sample_columns(
    file_path: Path,
    dataset: DatasetMetadata,
) -> None:
    """Select gene and sample columns, rename samples, and overwrite CSV."""
    if not file_path.is_file():
        raise FileNotFoundError(
            f"Dataset file does not exist: {file_path}"
        )

    df: pl.DataFrame = pl.read_csv(
        file_path,
        infer_schema_length=10_000,
    )

    if df.width == 0:
        raise ValueError(
            f"Dataset contains no columns: {file_path}"
        )

    gene_column: str = df.columns[0]

    sample_columns: list[str] = [
        *dataset.ctrl_cols,
        *dataset.ko_cols,
    ]

    missing_columns: list[str] = [
        column
        for column in sample_columns
        if column not in df.columns
    ]

    if missing_columns:
        raise ValueError(
            f"{dataset.dataset_name} is missing configured sample "
            f"columns: {missing_columns}\n"
            f"Available columns: {df.columns}"
        )

    duplicated_config_columns: set[str] = {
        column
        for column in sample_columns
        if sample_columns.count(column) > 1
    }

    if duplicated_config_columns:
        raise ValueError(
            f"{dataset.dataset_name} has duplicated sample names "
            f"in its TOML configuration: "
            f"{sorted(duplicated_config_columns)}"
        )

    control_expressions: list[pl.Expr] = [
        pl.col(original_name).alias(f"ctrl{index}")
        for index, original_name in enumerate(
            dataset.ctrl_cols,
            start=1,
        )
    ]

    knockout_expressions: list[pl.Expr] = [
        pl.col(original_name).alias(f"ko{index}")
        for index, original_name in enumerate(
            dataset.ko_cols,
            start=1,
        )
    ]

    standardized_df: pl.DataFrame = df.select(
        pl.col(gene_column).alias("gene"),
        *control_expressions,
        *knockout_expressions,
    )

    # Writing to a temporary file first avoids replacing the source
    # with a partially written file if an error occurs.
    temporary_file: Path = file_path.with_suffix(".tmp.csv")

    standardized_df.write_csv(temporary_file)
    temporary_file.replace(file_path)

    print(
        f"Saved {dataset.dataset_name}: "
        f"{standardized_df.height:,} rows, "
        f"{standardized_df.width} columns -> {file_path}"
    )


mouse_datasets: list[DatasetMetadata] = get_mouse_tec_low()

input_directory: Path = resolve_dir_from_root(
    "data/mtec_low/raw/intermediate/annotated_expression_matrix/mouse"
)

for dataset in mouse_datasets:
    # Keep this only if GSE65617 is processed separately.
    if dataset.dataset_name != "GSE2585":
        continue

    dataset_file: Path = (
        input_directory / f"{dataset.dataset_id}.csv"
    )

    standardize_sample_columns(
        file_path=dataset_file,
        dataset=dataset,
    )