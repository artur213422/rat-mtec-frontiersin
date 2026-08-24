from pathlib import Path
import re

import numpy as np
import polars as pl

from src.common.dataset_catalog.api.get_datasets import get_mtec_low
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.paths import resolve_dir_from_root


MetricRow = dict[str, str | int | float]


def get_expression_columns(
    data: pl.DataFrame,
    prefixes: tuple[str, ...],
) -> list[str]:
    """
    Finds columns such as:
    he1, he2, he3
    ctrl1, ctrl2, ctrl3
    control1, control2
    ko1, ko2, ko3
    """
    normalized_prefixes: set[str] = {
        prefix.casefold()
        for prefix in prefixes
    }

    pattern: re.Pattern[str] = re.compile(
        r"^(?P<prefix>[A-Za-z]+)_?(?P<index>\d+)$",
        flags=re.IGNORECASE,
    )

    matching_columns: list[tuple[int, str]] = []

    for column in data.columns:
        match: re.Match[str] | None = pattern.fullmatch(column)

        if match is None:
            continue

        prefix: str = match.group("prefix").casefold()

        if prefix not in normalized_prefixes:
            continue

        sample_index: int = int(match.group("index"))
        matching_columns.append((sample_index, column))

    matching_columns.sort(key=lambda item: item[0])

    return [
        column
        for _, column in matching_columns
    ]


def extract_numeric_matrix(
    data: pl.DataFrame,
    columns: list[str],
) -> np.ndarray:
    if len(columns) < 2:
        raise ValueError(
            f"At least two replicate columns are required. Found: {columns}"
        )

    matrix: np.ndarray = (
        data
        .select(
            [
                pl.col(column).cast(pl.Float64, strict=False)
                for column in columns
            ]
        )
        .to_numpy()
    )

    return np.asarray(matrix, dtype=np.float64)


def calculate_gene_standard_deviations(
    matrix: np.ndarray,
) -> np.ndarray:
    """
    Calculates one within-group SD for every gene.

    Missing values are allowed, but each gene must have at least
    two finite replicate values.
    """
    finite_mask: np.ndarray = np.isfinite(matrix)
    valid_counts: np.ndarray = finite_mask.sum(axis=1)

    row_sums: np.ndarray = np.where(
        finite_mask,
        matrix,
        0.0,
    ).sum(axis=1)

    row_means: np.ndarray = np.divide(
        row_sums,
        valid_counts,
        out=np.full(matrix.shape[0], np.nan),
        where=valid_counts > 0,
    )

    centered: np.ndarray = np.where(
        finite_mask,
        matrix - row_means[:, np.newaxis],
        0.0,
    )

    sum_squared_deviations: np.ndarray = np.square(
        centered
    ).sum(axis=1)

    gene_variances: np.ndarray = np.divide(
        sum_squared_deviations,
        valid_counts - 1,
        out=np.full(matrix.shape[0], np.nan),
        where=valid_counts >= 2,
    )

    gene_sds: np.ndarray = np.sqrt(gene_variances)

    return gene_sds[np.isfinite(gene_sds)]


def calculate_pairwise_correlations(
    matrix: np.ndarray,
) -> np.ndarray:
    correlations: list[float] = []

    sample_count: int = matrix.shape[1]

    for first_index in range(sample_count):
        for second_index in range(first_index + 1, sample_count):
            first_sample: np.ndarray = matrix[:, first_index]
            second_sample: np.ndarray = matrix[:, second_index]

            valid_mask: np.ndarray = (
                np.isfinite(first_sample)
                & np.isfinite(second_sample)
            )

            first_valid: np.ndarray = first_sample[valid_mask]
            second_valid: np.ndarray = second_sample[valid_mask]

            if first_valid.size < 2:
                continue

            if (
                np.std(first_valid) == 0.0
                or np.std(second_valid) == 0.0
            ):
                continue

            correlation: float = float(
                np.corrcoef(
                    first_valid,
                    second_valid,
                )[0, 1]
            )

            if np.isfinite(correlation):
                correlations.append(correlation)

    return np.asarray(correlations, dtype=np.float64)


def summarize_group_variation(
    dataset_id: str,
    group: str,
    matrix: np.ndarray,
) -> MetricRow:
    gene_sds: np.ndarray = calculate_gene_standard_deviations(
        matrix
    )

    correlations: np.ndarray = calculate_pairwise_correlations(
        matrix
    )

    if gene_sds.size == 0:
        raise ValueError(
            f"No valid gene SD values for {dataset_id}, group {group}."
        )

    mean_correlation: float = (
        float(np.mean(correlations))
        if correlations.size > 0
        else float("nan")
    )

    minimum_correlation: float = (
        float(np.min(correlations))
        if correlations.size > 0
        else float("nan")
    )

    return {
        "dataset_id": dataset_id,
        "group": group,
        "sample_count": int(matrix.shape[1]),
        "gene_count": int(gene_sds.size),
        "median_gene_sd": float(np.median(gene_sds)),
        "mean_gene_sd": float(np.mean(gene_sds)),
        "p90_gene_sd": float(np.percentile(gene_sds, 90)),
        "p95_gene_sd": float(np.percentile(gene_sds, 95)),
        "percent_sd_ge_0_5": float(
            100.0 * np.mean(gene_sds >= 0.5)
        ),
        "percent_sd_ge_1_0": float(
            100.0 * np.mean(gene_sds >= 1.0)
        ),
        "mean_pairwise_correlation": mean_correlation,
        "minimum_pairwise_correlation": minimum_correlation,
    }


def analyze_dataset_variation(
    dataset_id: str,
    data: pl.DataFrame,
) -> list[MetricRow]:
    control_columns: list[str] = get_expression_columns(
        data=data,
        prefixes=("he", "ctrl", "control"),
    )

    ko_columns: list[str] = get_expression_columns(
        data=data,
        prefixes=("ko",),
    )

    if not control_columns:
        raise ValueError(
            f"No HE/control columns found in {dataset_id}. "
            f"Available columns: {data.columns}"
        )

    if not ko_columns:
        raise ValueError(
            f"No KO columns found in {dataset_id}. "
            f"Available columns: {data.columns}"
        )

    print(f"{dataset_id}:")
    print(f"  HE/control columns: {control_columns}")
    print(f"  KO columns:         {ko_columns}")

    control_matrix: np.ndarray = extract_numeric_matrix(
        data=data,
        columns=control_columns,
    )

    ko_matrix: np.ndarray = extract_numeric_matrix(
        data=data,
        columns=ko_columns,
    )

    return [
        summarize_group_variation(
            dataset_id=dataset_id,
            group="HE",
            matrix=control_matrix,
        ),
        summarize_group_variation(
            dataset_id=dataset_id,
            group="KO",
            matrix=ko_matrix,
        ),
    ]


summary_file: Path = resolve_dir_from_root("data/mtec_low/analysis/general_statistics")/"dataset_variation.txt"


datasets: list[DatasetMetadata] = get_mtec_low()
variation_rows: list[MetricRow] = []

for dataset in datasets:
    print(f"\nRunning dataset: {dataset.dataset_id}")
    if dataset.dataset_id == "mouse_gse205410_mtec_low":
        continue

    tec_df: pl.DataFrame = dataset.get_standardized_df()

    dataset_rows: list[MetricRow] = analyze_dataset_variation(
        dataset_id=dataset.dataset_id,
        data=tec_df,
    )

    variation_rows.extend(dataset_rows)


variation_df: pl.DataFrame = (
    pl.DataFrame(variation_rows)
    .sort(
        ["dataset_id", "group"]
    )
)

print("\nDATASET VARIATION")
print(variation_df)

variation_df.write_csv(
    summary_file,
    separator="\t",
)

print(f"\nSaved variation summary to: {summary_file}")