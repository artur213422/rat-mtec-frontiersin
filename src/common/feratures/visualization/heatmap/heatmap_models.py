from dataclasses import dataclass
from typing import Self

import polars as pl

import numpy as np
import numpy.typing as npt

from src.common.types.enums import Colors


@dataclass(frozen=True, slots=True)
class HeatmapStyle:
    high_color: Colors
    low_color: Colors


    vcenter: float
    center_color: Colors | None = None
    nan_cell_color: Colors | None = None
    vmax: float | None = None
    vmin: float | None = None



@dataclass(slots=True)
class HeatmapData:
    row_labels: list[str]
    col_labels: list[str]
    matrix: npt.NDArray[np.float64]

    def move_col_to_left(self, col_label: str) -> None:
        if col_label not in self.col_labels:
            raise ValueError(f"Column label not found: {col_label}")

        old_index = self.col_labels.index(col_label)

        new_order = [old_index] + [
            i for i in range(len(self.col_labels))
            if i != old_index
        ]

        self.col_labels = [self.col_labels[i] for i in new_order]
        self.matrix = self.matrix[:, new_order]

    def sort_rows_by_col(self, col_label: str) -> None:
        if col_label not in self.col_labels:
            raise ValueError(f"Column label not found: {col_label}")

        col_index: int = self.col_labels.index(col_label)

        values: npt.NDArray[np.float64] = self.matrix[:, col_index]

        row_order: npt.NDArray[np.intp] = np.argsort(
            np.nan_to_num(values, nan=np.inf)
        )

        self.matrix = self.matrix[row_order, :]
        self.row_labels = [
            self.row_labels[i]
            for i in row_order
        ]

    @classmethod
    def from_long_frame(cls, *, data: pl.DataFrame, row_label_col: str, col_label_col: str, value_col: str,) -> Self:
        duplicate_pairs = (
            data
            .group_by(
                [row_label_col, col_label_col],
                maintain_order=True,
            )
            .len()
            .filter(pl.col("len") > 1)
        )
        if data.is_empty():
            raise ValueError("Cannot create heatmap data from empty DataFrame")

        if duplicate_pairs.height > 0:
            raise ValueError(f"Duplicate row/column pairs found: {duplicate_pairs}")

        pivoted = data.pivot(
            on=col_label_col,
            index=row_label_col,
            values=value_col,
            aggregate_function=None,
            maintain_order=True,
            sort_columns=False,
        )

        row_labels: list[str] = (
            pivoted
            .get_column(row_label_col)
            .to_list()
        )

        col_labels: list[str] = [
            col
            for col in pivoted.columns
            if col != row_label_col
        ]

        matrix: npt.NDArray[np.float64] = (
            pivoted
            .select(col_labels)
            .to_numpy()
            .astype(np.float64)
        )

        return cls(
            row_labels=row_labels,
            col_labels=col_labels,
            matrix=matrix,
        )


@dataclass(frozen=True, slots=True)
class HeatmapSpec:
    style: HeatmapStyle
    data: HeatmapData

