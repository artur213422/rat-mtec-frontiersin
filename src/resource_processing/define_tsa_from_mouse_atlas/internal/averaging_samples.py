import re

import polars as pl

def average_sample(matrix: pl.DataFrame, sample_columns: list[str], output_column_name: str) -> pl.DataFrame:
    if not sample_columns:
        raise ValueError(f"No replicate columns provided for group: {output_column_name}")

    if len(sample_columns) != len(set(sample_columns)):
        raise ValueError(f"Duplicate replicate columns provided for {output_column_name}: {sample_columns}")

    missing_columns: list[str] = [
        column
        for column in sample_columns
        if column not in matrix.columns
    ]

    if missing_columns:
        raise ValueError(f"Replicate columns missing for {output_column_name}: {missing_columns}")

    if output_column_name in matrix.columns and output_column_name not in sample_columns:
        raise ValueError(f"Output column already exists: {output_column_name}")

    non_numeric_columns: list[str] = [
        column
        for column in sample_columns
        if not matrix.schema[column].is_numeric()
    ]

    if non_numeric_columns:
        raise TypeError(f"Replicate columns must be numeric: {non_numeric_columns}")

    null_counts: dict[str, int] = (
        matrix
        .select(sample_columns)
        .null_count()
        .row(0, named=True)
    )

    columns_with_nulls: dict[str, int] = {
        column: count
        for column, count in null_counts.items()
        if count > 0
    }

    if columns_with_nulls:
        raise ValueError(f"Null values found in replicate columns for {output_column_name}: {columns_with_nulls}")

    sample_column_set: set[str] = set(sample_columns)

    # Place the averaged column where the first replicate occurred.
    first_sample_index: int = min(
        matrix.columns.index(column)
        for column in sample_columns
    )

    columns_before: list[str] = [
        column
        for column in matrix.columns[:first_sample_index]
        if column not in sample_column_set
    ]

    columns_after: list[str] = [
        column
        for column in matrix.columns[first_sample_index:]
        if column not in sample_column_set
           and column != output_column_name
    ]

    return matrix.select(
        *columns_before,
        pl.mean_horizontal(
            *[
                pl.col(column)
                for column in sample_columns
            ]
        ).alias(output_column_name),
        *columns_after,
    )