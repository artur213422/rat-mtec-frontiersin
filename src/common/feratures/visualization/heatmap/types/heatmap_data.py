from dataclasses import dataclass

import numpy as np
import polars as pl

import re

def natural_sort_key(value: str) -> list[str | int]:
    return [
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", value)
    ]

@dataclass(frozen=True, slots=True)
class HeatmapData:
    data: pl.DataFrame
    value_col: str
    col_col: str
    row_col: str

    def __post_init__(self) -> None:
        required_cols: set[str] = {
            self.row_col,
            self.col_col,
            self.value_col,
        }

        missing_cols: set[str] = required_cols - set(self.data.columns)

        if missing_cols:
            raise ValueError(f"Missing required columns: {sorted(missing_cols)}")

        if self.data.is_empty():
            raise ValueError("Cannot create heatmap from empty DataFrame.")

        row_labels: list[str] = (
            self.data
            .get_column(self.row_col)
            .unique()
            .cast(pl.String)
            .to_list()
        )

        ordered_rows = sorted(row_labels, key=natural_sort_key,)

        row_order: dict[str, int] = {
            label: index
            for index, label in enumerate(ordered_rows)
        }

        sorted_data = (
            self.data
            .with_columns(
                pl.col(self.row_col)
                .cast(pl.String)
                .replace_strict(row_order)
                .alias("__row_order")
            )
            .sort("__row_order")
            .drop("__row_order")
        )

        object.__setattr__(self, "data", sorted_data)

    def matrix_for(
        self,
        value_col: str,
        *,
        col_order: list[str] | None = None,
    ) -> np.ndarray[np.float64]:
        if value_col not in self.data.columns:
            raise ValueError(f"Column {value_col!r} not found in DataFrame.")

        pivoted = self.data.pivot(
            on=self.col_col,
            index=self.row_col,
            values=value_col,
            aggregate_function=None,
            maintain_order=True,
            sort_columns=False,
        )

        if col_order is not None:
            pivoted = pivoted.select(
                [
                    self.row_col,
                    *col_order,
                ]
            )

        return (
            pivoted
            .drop(self.row_col)
            .to_numpy()
            .astype(np.float64)
        )

    @property
    def row_labels(self) -> list[str]:
        return (
            self.data
            .get_column(self.row_col)
            .unique(maintain_order=True)
            .cast(pl.String)
            .to_list()
        )

    @property
    def col_labels(self) -> list[str]:
        return (
            self.data
            .get_column(self.col_col)
            .unique(maintain_order=True)
            .cast(pl.String)
            .to_list()
        )