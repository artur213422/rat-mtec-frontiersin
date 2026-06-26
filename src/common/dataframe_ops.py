import re
from pathlib import Path

import polars as pl

from src.common.models.dataset_metadata.model import DatasetMetadata
from src.common.types.enums import Genotypes

ID_COLUMNS: tuple[str, str] = ("gene_symbol", "ortholog_id")


def get_genotype_expression_values(row: dict[str, object], gt: Genotypes) -> list[float]:
    pattern = re.compile(rf"^{gt.value}\d+$")

    return [
        float(value)
        for col, value in row.items()
        if pattern.fullmatch(col) and value is not None
    ]


def _get_genotype_columns(df: pl.DataFrame, gt: Genotypes) -> list[str]:
    pattern: re.Pattern[str] = re.compile(rf"^{re.escape(gt.value)}\d+$")

    return [
        column
        for column in df.columns
        if pattern.fullmatch(column) is not None
    ]

def _build_dataset_frame(metadata: DatasetMetadata, ctrl_start_index: int, ko_start_index: int) -> tuple[pl.DataFrame, int, int, list[str], list[str]]:
    file: Path = metadata.get_standardized_data_file()
    df: pl.DataFrame = pl.read_csv(file)

    missing_id_columns: list[str] = [
        column for column in ID_COLUMNS
        if column not in df.columns
    ]

    if missing_id_columns:
        raise ValueError(
            f"{metadata.dataset_id}: missing columns: {missing_id_columns}"
        )

    ctrl_cols: list[str] = _get_genotype_columns(df, Genotypes.CTRL)
    ko_cols: list[str] = _get_genotype_columns(df, Genotypes.KO)

    if not ctrl_cols:
        raise ValueError(f"{metadata.dataset_id}: no CTRL columns found")

    if not ko_cols:
        raise ValueError(f"{metadata.dataset_id}: no KO columns found")

    rename_map: dict[str, str] = {}

    new_ctrl_cols: list[str] = []
    for old_col in ctrl_cols:
        new_col: str = f"ctrl{ctrl_start_index}"
        rename_map[old_col] = new_col
        new_ctrl_cols.append(new_col)
        ctrl_start_index += 1

    new_ko_cols: list[str] = []
    for old_col in ko_cols:
        new_col = f"ko{ko_start_index}"
        rename_map[old_col] = new_col
        new_ko_cols.append(new_col)
        ko_start_index += 1

    selected_df: pl.DataFrame = (
        df
        .select([*ID_COLUMNS, *ctrl_cols, *ko_cols])
        .rename(rename_map)
    )

    return selected_df, ctrl_start_index, ko_start_index, new_ctrl_cols, new_ko_cols


def combine_datasets_with_nulls(datasets: list[DatasetMetadata]) -> pl.DataFrame:
    if not datasets:
        return pl.DataFrame(
            schema={
                "gene_symbol": pl.Utf8,
                "ortholog_id": pl.Utf8,
            }
        )

    frames: list[pl.DataFrame] = []

    ctrl_index: int = 1
    ko_index: int = 1

    all_ctrl_cols: list[str] = []
    all_ko_cols: list[str] = []

    for metadata in datasets:
        frame, ctrl_index, ko_index, ctrl_cols, ko_cols = _build_dataset_frame(
            metadata=metadata,
            ctrl_start_index=ctrl_index,
            ko_start_index=ko_index,
        )

        frames.append(frame)
        all_ctrl_cols.extend(ctrl_cols)
        all_ko_cols.extend(ko_cols)

    combined: pl.DataFrame = frames[0]

    for frame in frames[1:]:
        combined = combined.join(
            frame,
            on=list(ID_COLUMNS),
            how="full",
            coalesce=True,
        )

    print_null_gene_summary(combined)
    return combined.select([*ID_COLUMNS, *all_ctrl_cols, *all_ko_cols])

def print_null_gene_summary(data: pl.DataFrame) -> None:
    sample_cols: list[str] = [
        column
        for column in data.columns
        if column not in ID_COLUMNS
    ]

    non_null_value_count: pl.Expr = pl.sum_horizontal(
        pl.col(sample_cols).is_not_null().cast(pl.Int64)
    )

    genes_with_nulls: int = data.filter(
        pl.any_horizontal(pl.col(sample_cols).is_null())
    ).height

    fully_populated_genes: int = data.filter(
        pl.all_horizontal(pl.col(sample_cols).is_not_null())
    ).height

    genes_with_less_than_5_values: int = data.filter(
        non_null_value_count < 5
    ).height

    total_genes: int = data.height

    print(f"Total genes: {total_genes}")
    print(f"Genes containing nulls: {genes_with_nulls}")
    print(f"Fully populated genes: {fully_populated_genes}")
    print(f"Genes with less than 5 values: {genes_with_less_than_5_values}")