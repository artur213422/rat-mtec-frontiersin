import math
import statistics
from typing import Sequence

import polars as pl

from src.common.dataframe_ops import get_genotype_expression_values
from src.common.types.enums import Genotypes




def calc_welch_t_rank(control_values: Sequence[float], case_values: Sequence[float]) -> float | None:
    if len(control_values) < 2:
        return None

    if len(case_values) < 2:
        return None

    control_mean = statistics.fmean(control_values)
    case_mean = statistics.fmean(case_values)

    control_variance = statistics.variance(control_values)
    case_variance = statistics.variance(case_values)

    standard_error = math.sqrt(
        (control_variance / len(control_values))
        + (case_variance / len(case_values))
    )

    if standard_error == 0:
        return None

    return (case_mean - control_mean) / standard_error


def calculate_rnk_matrix(df: pl.DataFrame) -> pl.DataFrame:
    if "gene_symbol" not in df.columns:
        raise ValueError("Missing required column: gene_symbol")

    rows: list[dict[str, str | float]] = []

    for row in df.iter_rows(named=True):
        gene_symbol = row.get("gene_symbol")

        if not isinstance(gene_symbol, str):
            continue

        control_values = get_genotype_expression_values(row, Genotypes.KO)
        case_values = get_genotype_expression_values(row, Genotypes.CTRL)

        rank = calc_welch_t_rank(
            control_values=control_values,
            case_values=case_values,
        )

        if rank is None:
            continue

        if not math.isfinite(rank):
            continue

        rows.append(
            {
                "gene_symbol": gene_symbol,
                "rank": rank,
            }
        )

    return (
        pl.DataFrame(
            rows,
            schema={
                "gene_symbol": pl.String,
                "rank": pl.Float64,
            },
        )
        .unique(subset=["gene_symbol"], keep="first")
        .sort("rank", descending=True)
    )



def calculate_rnk_matrix_rat_special(df: pl.DataFrame) -> pl.DataFrame:
    required_columns: set[str] = {
        "has_clean_ortholog",
        "gene_symbol",
        "mouse_ortholog_name",
    }

    missing_columns: set[str] = required_columns - set(df.columns)

    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    df = (
        df
        .filter(pl.col("has_clean_ortholog"))
        .drop("gene_symbol")
        .rename({"mouse_ortholog_name": "gene_symbol"})
    )

    rows: list[dict[str, str | float]] = []

    for row in df.iter_rows(named=True):
        gene_symbol = row.get("gene_symbol")

        if not isinstance(gene_symbol, str):
            continue

        control_values = get_genotype_expression_values(
            row,
            Genotypes.CTRL,
        )
        case_values = get_genotype_expression_values(
            row,
            Genotypes.KO,
        )

        rank = calc_welch_t_rank(
            control_values=control_values,
            case_values=case_values,
        )

        if rank is None:
            continue

        if not math.isfinite(rank):
            continue

        rows.append(
            {
                "gene_symbol": gene_symbol,
                "rank": rank,
            }
        )

    return (
        pl.DataFrame(
            rows,
            schema={
                "gene_symbol": pl.String,
                "rank": pl.Float64,
            },
        )
        .unique(subset=["gene_symbol"], keep="first")
        .sort("rank", descending=True)
    )
