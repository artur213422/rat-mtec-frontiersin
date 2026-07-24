from __future__ import annotations

import math
from statistics import fmean, mean

import polars as pl

from src.common.dataframe_ops import get_genotype_expression_values
from src.common.types.enums import Genotypes

def calculate_log2_fold_change(ctrl_values: list[float], ko_values: list[float]) -> float:
    if not ctrl_values:
        raise ValueError("Control values are empty.")

    if not ko_values:
        raise ValueError("KO values are empty.")

    return mean(ko_values) - mean(ctrl_values)

def calculate_log2_mean(values: list[float]) -> float | None:
    if not values:
        return None

    return math.log2(fmean(values) + 1.0)


def calculate_fold_changes(data: pl.DataFrame) -> pl.DataFrame:
    if data.is_empty():
        raise ValueError("Cannot calculate fold changes for an empty DataFrame.")

    ctrl_log2_means: list[float] = []
    ko_log2_means: list[float] = []
    log2_fold_changes: list[float] = []

    for row_number, row in enumerate(data.iter_rows(named=True), start=1):
        controls: list[float] = get_genotype_expression_values(
            row,
            Genotypes.CTRL,
        )
        knockouts: list[float] = get_genotype_expression_values(
            row,
            Genotypes.KO,
        )

        if not controls:
            raise ValueError(f"No CTRL expression values in row {row_number}.")

        if not knockouts:
            raise ValueError(f"No KO expression values in row {row_number}.")

        ctrl_log2_mean: float = fmean(controls)
        ko_log2_mean: float = fmean(knockouts)

        if not math.isfinite(ctrl_log2_mean):
            raise ValueError( f"Non-finite CTRL mean in row {row_number}.")

        if not math.isfinite(ko_log2_mean):
            raise ValueError(
                f"Non-finite KO mean in row {row_number}."
            )

        log2_fold_change: float = (
            ko_log2_mean - ctrl_log2_mean
        )

        ctrl_log2_means.append(ctrl_log2_mean)
        ko_log2_means.append(ko_log2_mean)
        log2_fold_changes.append(log2_fold_change)

    return data.with_columns(
        pl.Series(
            name="ctrl_log2_mean",
            values=ctrl_log2_means,
            dtype=pl.Float64,
        ),
        pl.Series(
            name="ko_log2_mean",
            values=ko_log2_means,
            dtype=pl.Float64,
        ),
        pl.Series(
            name="log2_fold_change",
            values=log2_fold_changes,
            dtype=pl.Float64,
        ),
    )