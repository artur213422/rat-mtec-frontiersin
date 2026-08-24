from __future__ import annotations

from itertools import combinations
from pathlib import Path

import numpy as np
import numpy.typing as npt
import polars as pl
from scipy.stats import (
    chi2,
    kendalltau,
    rankdata,
    spearmanr,
)

from src.common.paths import (
    resolve_dir_from_root,
    try_find_file_in_project,
)


##########################################################################################
#                                       PREPARE DATA
##########################################################################################

TOP_GENE_COUNT: int = 1_000

analysis_folder: Path = resolve_dir_from_root(
    "data/mtec_low/analysis/affy_vs_rna_seq"
)

analysis_folder.mkdir(
    parents=True,
    exist_ok=True,
)

output_file: Path = (
    analysis_folder
    / f"summary_top_{TOP_GENE_COUNT}_consensus_effects.txt"
)

analysed_file: Path = try_find_file_in_project(
    "union_genes.csv"
)

EXPECTED_DATASETS: set[str] = {
    "mouse_gse2585_mtec_low",
    "mouse_gse65617_mtec_low",
    "mouse_gse205410_mtec_low",
}

REQUIRED_COLUMNS: set[str] = {
    "gene_symbol",
    "log2_fold_change",
    "dataset_id",
}


##########################################################################################
#                                       HELPERS
##########################################################################################


def format_number(
    value: float,
    *,
    decimal_places: int = 4,
) -> str:
    if not np.isfinite(value):
        return "NA"

    return f"{value:.{decimal_places}f}"


def format_p_value(value: float) -> str:
    if not np.isfinite(value):
        return "NA"

    if value == 0.0:
        return "<1e-300"

    if value < 0.0001:
        return f"{value:.3e}"

    return f"{value:.4f}"


def format_percentage(
    numerator: int,
    denominator: int,
    *,
    decimal_places: int = 1,
) -> str:
    if denominator <= 0:
        return "NA"

    percentage: float = (
        100.0 * numerator / denominator
    )

    return f"{percentage:.{decimal_places}f}%"


def calculate_kendalls_w(
    values: npt.NDArray[np.float64],
) -> tuple[float, float, float]:
    """
    Calculate Kendall's coefficient of concordance.

    Parameters
    ----------
    values
        Matrix with genes in rows and datasets in columns.

    Returns
    -------
    tuple
        Kendall's W, chi-square statistic, and approximate p-value.

    Notes
    -----
    Ranks are assigned independently within each dataset.

    Tied values receive average ranks. A correction for tied groups is
    included in the denominator.
    """
    if values.ndim != 2:
        raise ValueError(
            "Kendall's W input must be a two-dimensional matrix."
        )

    n_genes, n_datasets = values.shape

    if n_genes < 2:
        raise ValueError(
            "At least two genes are required for Kendall's W."
        )

    if n_datasets < 2:
        raise ValueError(
            "At least two datasets are required for Kendall's W."
        )

    ranked_values: npt.NDArray[np.float64] = np.column_stack(
        [
            rankdata(
                values[:, dataset_index],
                method="average",
            )
            for dataset_index in range(n_datasets)
        ]
    ).astype(np.float64)

    rank_sums: npt.NDArray[np.float64] = (
        ranked_values.sum(axis=1)
    )

    expected_rank_sum: float = (
        n_datasets * (n_genes + 1) / 2.0
    )

    squared_deviation_sum: float = float(
        np.square(
            rank_sums - expected_rank_sum
        ).sum()
    )

    total_tie_correction: float = 0.0

    for dataset_index in range(n_datasets):
        _, tie_counts = np.unique(
            values[:, dataset_index],
            return_counts=True,
        )

        tied_groups: npt.NDArray[np.int64] = (
            tie_counts[tie_counts > 1]
        )

        total_tie_correction += float(
            np.sum(
                tied_groups**3 - tied_groups
            )
        )

    denominator: float = (
        n_datasets**2
        * (n_genes**3 - n_genes)
        - n_datasets * total_tie_correction
    )

    if denominator <= 0:
        raise ValueError(
            "Kendall's W denominator is not positive. "
            "One or more datasets may contain constant values."
        )

    kendalls_w: float = (
        12.0 * squared_deviation_sum
        / denominator
    )

    kendalls_w = float(
        np.clip(
            kendalls_w,
            0.0,
            1.0,
        )
    )

    chi_square_statistic: float = (
        n_datasets
        * (n_genes - 1)
        * kendalls_w
    )

    p_value: float = float(
        chi2.sf(
            chi_square_statistic,
            df=n_genes - 1,
        )
    )

    return (
        kendalls_w,
        chi_square_statistic,
        p_value,
    )


def get_complete_shared_genes(
    data: pl.DataFrame,
    *,
    dataset_ids: list[str],
) -> pl.DataFrame:
    """
    Retain genes with one finite fold-change value in every dataset.
    """
    number_of_datasets: int = len(dataset_ids)

    valid_data: pl.DataFrame = data.filter(
        pl.col("log2_fold_change").is_not_null()
        & pl.col("log2_fold_change").is_not_nan()
        & pl.col("log2_fold_change").is_finite()
    )

    shared_gene_symbols: pl.DataFrame = (
        valid_data
        .group_by("gene_symbol")
        .agg(
            pl.col("dataset_id")
            .n_unique()
            .alias("n_datasets")
        )
        .filter(
            pl.col("n_datasets")
            == number_of_datasets
        )
        .select("gene_symbol")
    )

    return (
        valid_data
        .join(
            shared_gene_symbols,
            on="gene_symbol",
            how="semi",
        )
        .sort(
            "gene_symbol",
            "dataset_id",
        )
    )


def select_consensus_top_effect_genes(
    wide_frame: pl.DataFrame,
    *,
    dataset_ids: list[str],
    top_n: int,
) -> pl.DataFrame:
    """
    Select genes with the strongest relative effects across datasets.

    Absolute log2 fold changes are ranked independently within each
    dataset. Rank 1 is the strongest absolute effect in that dataset.

    The mean absolute-effect rank is then calculated across datasets.
    Genes with the lowest mean rank are retained.

    This avoids penalizing Affymetrix for having a smaller fold-change
    range than RNA-seq.
    """
    if top_n <= 0:
        raise ValueError(
            "top_n must be positive."
        )

    if top_n > wide_frame.height:
        raise ValueError(
            f"Requested {top_n:,} genes, but only "
            f"{wide_frame.height:,} shared genes are available."
        )

    rank_columns: list[str] = []
    rank_expressions: list[pl.Expr] = []

    for dataset_id in dataset_ids:
        rank_column: str = (
            f"{dataset_id}__absolute_effect_rank"
        )

        rank_columns.append(rank_column)

        rank_expressions.append(
            pl.col(dataset_id)
            .abs()
            .rank(
                method="average",
                descending=True,
            )
            .alias(rank_column)
        )

    ranked_frame: pl.DataFrame = (
        wide_frame
        .with_columns(rank_expressions)
        .with_columns(
            pl.mean_horizontal(
                *[
                    pl.col(rank_column)
                    for rank_column in rank_columns
                ]
            ).alias(
                "mean_absolute_effect_rank"
            ),

            pl.max_horizontal(
                *[
                    pl.col(rank_column)
                    for rank_column in rank_columns
                ]
            ).alias(
                "worst_absolute_effect_rank"
            ),
        )
    )

    return (
        ranked_frame
        .sort(
            [
                "mean_absolute_effect_rank",
                "worst_absolute_effect_rank",
                "gene_symbol",
            ]
        )
        .head(top_n)
    )


##########################################################################################
#                                       READ DATA
##########################################################################################

source_df: pl.DataFrame = pl.read_csv(
    analysed_file,
    infer_schema_length=10_000,
    schema_overrides={
        "gene_symbol": pl.String,
        "dataset_id": pl.String,
        "log2_fold_change": pl.Float64,
    },
    null_values=[
        "",
        "NA",
        "NaN",
        "nan",
        "NULL",
        "null",
    ],
)

missing_columns: set[str] = (
    REQUIRED_COLUMNS - set(source_df.columns)
)

if missing_columns:
    raise ValueError(
        f"Input file is missing required columns: "
        f"{sorted(missing_columns)}"
    )

source_df = (
    source_df
    .select(
        pl.col("gene_symbol")
        .cast(pl.String)
        .str.strip_chars(),

        pl.col("log2_fold_change")
        .cast(pl.Float64),

        pl.col("dataset_id")
        .cast(pl.String)
        .str.strip_chars(),
    )
    .filter(
        pl.col("gene_symbol").is_not_null()
        & (pl.col("gene_symbol") != "")
        & pl.col("dataset_id").is_not_null()
        & (pl.col("dataset_id") != "")
    )
)


##########################################################################################
#                                   VALIDATE DATASETS
##########################################################################################

observed_datasets: set[str] = set(
    source_df
    .get_column("dataset_id")
    .unique()
    .to_list()
)

if observed_datasets != EXPECTED_DATASETS:
    missing_datasets: set[str] = (
        EXPECTED_DATASETS - observed_datasets
    )

    unexpected_datasets: set[str] = (
        observed_datasets - EXPECTED_DATASETS
    )

    raise ValueError(
        "Dataset IDs do not match the expected datasets.\n"
        f"Missing: {sorted(missing_datasets)}\n"
        f"Unexpected: {sorted(unexpected_datasets)}"
    )

dataset_ids: list[str] = sorted(
    observed_datasets
)

duplicate_pairs_df: pl.DataFrame = (
    source_df
    .group_by(
        "gene_symbol",
        "dataset_id",
    )
    .len()
    .filter(
        pl.col("len") > 1
    )
)

if duplicate_pairs_df.height > 0:
    raise ValueError(
        "Duplicate gene/dataset combinations were found:\n"
        f"{duplicate_pairs_df.head(20)}"
    )

dataset_gene_counts: dict[str, int] = {
    str(row["dataset_id"]): int(row["n_genes"])
    for row in (
        source_df
        .group_by("dataset_id")
        .agg(
            pl.col("gene_symbol")
            .n_unique()
            .alias("n_genes")
        )
        .sort("dataset_id")
        .iter_rows(named=True)
    )
}

total_unique_genes: int = (
    source_df
    .get_column("gene_symbol")
    .n_unique()
)


##########################################################################################
#                                 EXTRACT COMPLETE OVERLAP
##########################################################################################

shared_df: pl.DataFrame = get_complete_shared_genes(
    source_df,
    dataset_ids=dataset_ids,
)

shared_gene_count: int = (
    shared_df
    .get_column("gene_symbol")
    .n_unique()
)

expected_shared_rows: int = (
    shared_gene_count
    * len(dataset_ids)
)

if shared_df.height != expected_shared_rows:
    raise ValueError(
        f"Expected {expected_shared_rows:,} shared rows, "
        f"but found {shared_df.height:,}."
    )

wide_df: pl.DataFrame = (
    shared_df
    .pivot(
        on="dataset_id",
        index="gene_symbol",
        values="log2_fold_change",
        aggregate_function=None,
        maintain_order=True,
        sort_columns=False,
    )
    .select(
        "gene_symbol",
        *dataset_ids,
    )
    .sort("gene_symbol")
)

complete_effect_matrix: npt.NDArray[np.float64] = (
    wide_df
    .select(dataset_ids)
    .to_numpy()
    .astype(np.float64)
)

if not np.isfinite(
    complete_effect_matrix
).all():
    raise ValueError(
        "The shared-gene matrix contains non-finite values."
    )


##########################################################################################
#                             SELECT TOP CONSENSUS EFFECTS
##########################################################################################

top_effects_df: pl.DataFrame = (
    select_consensus_top_effect_genes(
        wide_df,
        dataset_ids=dataset_ids,
        top_n=TOP_GENE_COUNT,
    )
)

selected_gene_count: int = (
    top_effects_df.height
)

if selected_gene_count != TOP_GENE_COUNT:
    raise ValueError(
        f"Expected {TOP_GENE_COUNT:,} selected genes, "
        f"but found {selected_gene_count:,}."
    )

effect_matrix: npt.NDArray[np.float64] = (
    top_effects_df
    .select(dataset_ids)
    .to_numpy()
    .astype(np.float64)
)

if not np.isfinite(effect_matrix).all():
    raise ValueError(
        "The selected top-effect matrix contains "
        "non-finite values."
    )

for dataset_index, dataset_id in enumerate(
    dataset_ids
):
    number_of_unique_values: int = np.unique(
        effect_matrix[:, dataset_index]
    ).size

    if number_of_unique_values < 2:
        raise ValueError(
            f"Dataset '{dataset_id}' contains no "
            "fold-change variation among selected genes."
        )


##########################################################################################
#                               SELECTION DIAGNOSTICS
##########################################################################################

individual_top_counts: dict[str, int] = {}

for dataset_id in dataset_ids:
    rank_column: str = (
        f"{dataset_id}__absolute_effect_rank"
    )

    individual_top_counts[dataset_id] = (
        top_effects_df
        .filter(
            pl.col(rank_column)
            <= TOP_GENE_COUNT
        )
        .height
    )

all_datasets_top_count: int = (
    top_effects_df
    .filter(
        pl.all_horizontal(
            *[
                pl.col(
                    f"{dataset_id}"
                    "__absolute_effect_rank"
                )
                <= TOP_GENE_COUNT
                for dataset_id in dataset_ids
            ]
        )
    )
    .height
)


##########################################################################################
#                               PAIRWISE RANK ANALYSIS
##########################################################################################

pairwise_results: list[
    dict[str, str | float | int]
] = []

for first_dataset, second_dataset in combinations(
    dataset_ids,
    2,
):
    first_values: npt.NDArray[np.float64] = (
        top_effects_df
        .get_column(first_dataset)
        .to_numpy()
        .astype(np.float64)
    )

    second_values: npt.NDArray[np.float64] = (
        top_effects_df
        .get_column(second_dataset)
        .to_numpy()
        .astype(np.float64)
    )

    spearman_result = spearmanr(
        first_values,
        second_values,
        nan_policy="raise",
    )

    kendall_result = kendalltau(
        first_values,
        second_values,
        nan_policy="raise",
        variant="b",
    )

    nonzero_mask: npt.NDArray[np.bool_] = (
        (first_values != 0.0)
        & (second_values != 0.0)
    )

    nonzero_pair_count: int = int(
        nonzero_mask.sum()
    )

    first_signs: npt.NDArray[np.float64] = np.sign(
        first_values[nonzero_mask]
    )

    second_signs: npt.NDArray[np.float64] = np.sign(
        second_values[nonzero_mask]
    )

    same_direction_count: int = int(
        (
            first_signs == second_signs
        ).sum()
    )

    opposite_direction_count: int = int(
        (
            first_signs != second_signs
        ).sum()
    )

    pairwise_results.append(
        {
            "first_dataset": first_dataset,
            "second_dataset": second_dataset,
            "spearman_rho": float(
                spearman_result.statistic
            ),
            "spearman_p_value": float(
                spearman_result.pvalue
            ),
            "kendall_tau": float(
                kendall_result.statistic
            ),
            "kendall_p_value": float(
                kendall_result.pvalue
            ),
            "nonzero_pair_count": nonzero_pair_count,
            "same_direction_count": same_direction_count,
            "opposite_direction_count": (
                opposite_direction_count
            ),
        }
    )


##########################################################################################
#                              THREE-DATASET DIRECTION
##########################################################################################

all_negative_mask: npt.NDArray[np.bool_] = np.all(
    effect_matrix < 0.0,
    axis=1,
)

all_positive_mask: npt.NDArray[np.bool_] = np.all(
    effect_matrix > 0.0,
    axis=1,
)

mixed_direction_mask: npt.NDArray[np.bool_] = (
    np.any(
        effect_matrix < 0.0,
        axis=1,
    )
    & np.any(
        effect_matrix > 0.0,
        axis=1,
    )
)

contains_zero_mask: npt.NDArray[np.bool_] = np.any(
    effect_matrix == 0.0,
    axis=1,
)

unanimous_down_count: int = int(
    all_negative_mask.sum()
)

unanimous_up_count: int = int(
    all_positive_mask.sum()
)

mixed_direction_count: int = int(
    mixed_direction_mask.sum()
)

contains_zero_count: int = int(
    contains_zero_mask.sum()
)

unanimous_direction_count: int = (
    unanimous_down_count
    + unanimous_up_count
)


##########################################################################################
#                              GLOBAL RANK CONCORDANCE
##########################################################################################

(
    kendalls_w,
    kendalls_w_chi_square,
    kendalls_w_p_value,
) = calculate_kendalls_w(
    effect_matrix
)

spearman_values: list[float] = [
    float(result["spearman_rho"])
    for result in pairwise_results
]

kendall_values: list[float] = [
    float(result["kendall_tau"])
    for result in pairwise_results
]

unrelated_w_baseline: float = (
    1.0 / len(dataset_ids)
)


##########################################################################################
#                                       BUILD REPORT
##########################################################################################

report_lines: list[str] = [
    "MTEC-LOW TOP-EFFECT RANK-CONCORDANCE REPORT",
    "=" * 76,
    "",
    "Analysis objective",
    "-" * 76,
    (
        "Compare the direction and ordering of the strongest "
        "KO-versus-control effects across datasets without requiring "
        "equal fold-change magnitudes."
    ),
    (
        f"The analysis retains {TOP_GENE_COUNT:,} genes selected using "
        "consensus absolute-effect ranks."
    ),
    (
        "Absolute log2 fold changes were ranked independently within each "
        "dataset. Rank 1 represents the strongest effect in that dataset."
    ),
    (
        "The mean absolute-effect rank across all three datasets was used "
        "to select genes."
    ),
    "",
    "Input",
    "-" * 76,
    f"Source file: {analysed_file}",
    f"Input rows: {source_df.height:,}",
    f"Unique genes in input: {total_unique_genes:,}",
    f"Datasets: {len(dataset_ids)}",
]

for dataset_id in dataset_ids:
    report_lines.append(
        f"  {dataset_id}: "
        f"{dataset_gene_counts[dataset_id]:,} genes"
    )

report_lines.extend(
    [
        "",
        "Complete shared-gene set",
        "-" * 76,
        (
            "Genes with one finite effect in every dataset: "
            f"{shared_gene_count:,}"
        ),
        (
            "Genes excluded because they were incomplete: "
            f"{total_unique_genes - shared_gene_count:,}"
        ),
        "",
        "Strong-effect gene selection",
        "-" * 76,
        (
            f"Genes retained for analysis: "
            f"{selected_gene_count:,}"
        ),
        (
            "Selection method: mean within-dataset rank of absolute "
            "log2 fold change."
        ),
        (
            "A gene can therefore be selected whether it is strongly "
            "decreased or strongly increased."
        ),
        "",
        "Selected genes individually ranked in each dataset's top set",
    ]
)

for dataset_id in dataset_ids:
    count: int = individual_top_counts[dataset_id]

    report_lines.append(
        f"  {dataset_id}: "
        f"{count:,}/{selected_gene_count:,} "
        f"({format_percentage(count, selected_gene_count)})"
    )

report_lines.extend(
    [
        (
            "Selected genes ranked among the strongest "
            f"{TOP_GENE_COUNT:,} in all three datasets: "
            f"{all_datasets_top_count:,}/"
            f"{selected_gene_count:,} "
            f"({format_percentage(all_datasets_top_count, selected_gene_count)})"
        ),
        "",
        "Pairwise rank and direction agreement",
        "-" * 76,
    ]
)

for result in pairwise_results:
    first_dataset: str = str(
        result["first_dataset"]
    )

    second_dataset: str = str(
        result["second_dataset"]
    )

    nonzero_pair_count = int(
        result["nonzero_pair_count"]
    )

    same_direction_count = int(
        result["same_direction_count"]
    )

    opposite_direction_count = int(
        result["opposite_direction_count"]
    )

    report_lines.extend(
        [
            first_dataset,
            f"versus {second_dataset}",
            (
                "  Spearman rho: "
                f"{format_number(float(result['spearman_rho']))}"
            ),
            (
                "  Spearman p-value: "
                f"{format_p_value(float(result['spearman_p_value']))}"
            ),
            (
                "  Kendall tau-b: "
                f"{format_number(float(result['kendall_tau']))}"
            ),
            (
                "  Kendall p-value: "
                f"{format_p_value(float(result['kendall_p_value']))}"
            ),
            (
                "  Genes with non-zero effects in both: "
                f"{nonzero_pair_count:,}"
            ),
            (
                "  Same fold-change direction: "
                f"{same_direction_count:,}/"
                f"{nonzero_pair_count:,} "
                f"({format_percentage(same_direction_count, nonzero_pair_count)})"
            ),
            (
                "  Opposite fold-change direction: "
                f"{opposite_direction_count:,}/"
                f"{nonzero_pair_count:,} "
                f"({format_percentage(opposite_direction_count, nonzero_pair_count)})"
            ),
            "",
        ]
    )

report_lines.extend(
    [
        "Pairwise correlation summary",
        "-" * 76,
        (
            "Mean Spearman rho: "
            f"{format_number(float(np.mean(spearman_values)))}"
        ),
        (
            "Median Spearman rho: "
            f"{format_number(float(np.median(spearman_values)))}"
        ),
        (
            "Spearman rho range: "
            f"{format_number(float(np.min(spearman_values)))}"
            " to "
            f"{format_number(float(np.max(spearman_values)))}"
        ),
        (
            "Mean Kendall tau-b: "
            f"{format_number(float(np.mean(kendall_values)))}"
        ),
        (
            "Median Kendall tau-b: "
            f"{format_number(float(np.median(kendall_values)))}"
        ),
        "",
        "Three-dataset directional agreement",
        "-" * 76,
        (
            "Decreased in all three datasets: "
            f"{unanimous_down_count:,}/"
            f"{selected_gene_count:,} "
            f"({format_percentage(unanimous_down_count, selected_gene_count)})"
        ),
        (
            "Increased in all three datasets: "
            f"{unanimous_up_count:,}/"
            f"{selected_gene_count:,} "
            f"({format_percentage(unanimous_up_count, selected_gene_count)})"
        ),
        (
            "Same direction in all three datasets: "
            f"{unanimous_direction_count:,}/"
            f"{selected_gene_count:,} "
            f"({format_percentage(unanimous_direction_count, selected_gene_count)})"
        ),
        (
            "Contains both positive and negative directions: "
            f"{mixed_direction_count:,}/"
            f"{selected_gene_count:,} "
            f"({format_percentage(mixed_direction_count, selected_gene_count)})"
        ),
        (
            "Contains at least one exact zero effect: "
            f"{contains_zero_count:,}/"
            f"{selected_gene_count:,} "
            f"({format_percentage(contains_zero_count, selected_gene_count)})"
        ),
        "",
        "Global three-dataset rank concordance",
        "-" * 76,
        (
            "Kendall's W: "
            f"{format_number(kendalls_w)}"
        ),
        (
            "Approximate unrelated-ranking baseline for "
            f"{len(dataset_ids)} datasets: "
            f"{format_number(unrelated_w_baseline)}"
        ),
        (
            "Approximate chi-square statistic: "
            f"{format_number(kendalls_w_chi_square)}"
        ),
        (
            "Approximate p-value: "
            f"{format_p_value(kendalls_w_p_value)}"
        ),
        "",
        "Interpretation notes",
        "-" * 76,
        (
            "Spearman rho measures whether two datasets order the selected "
            "genes similarly from the most negative to the most positive "
            "fold change."
        ),
        (
            "Kendall tau-b is another pairwise rank-agreement measure and "
            "accounts for tied values."
        ),
        (
            "Same-direction percentages ask only whether fold-change signs "
            "agree; they do not require equal fold-change magnitudes."
        ),
        (
            "Kendall's W summarizes common rank agreement across all three "
            "datasets. A value of 1 indicates identical rankings."
        ),
        (
            f"With {len(dataset_ids)} datasets, unrelated rankings have an "
            f"expected W near {unrelated_w_baseline:.3f}; W should therefore "
            "not be interpreted as a direct percentage."
        ),
        (
            "Genes were selected from the same effect values later used for "
            "correlation analysis. Correlation p-values in this report are "
            "therefore descriptive after selection, not independent formal "
            "tests."
        ),
        (
            "Adjusted differential-expression p-values were not used for "
            "selection because statistical power differs substantially "
            "between datasets."
        ),
        (
            "A near-zero fold change means little KO-versus-control change. "
            "It does not necessarily mean that the gene has low absolute "
            "expression."
        ),
        "",
    ]
)


##########################################################################################
#                                       SAVE REPORT
##########################################################################################

output_file.write_text(
    "\n".join(report_lines),
    encoding="utf-8",
)

print(f"Saved: {output_file}")
print(
    f"Complete shared genes: "
    f"{shared_gene_count:,}"
)
print(
    f"Top consensus-effect genes analyzed: "
    f"{selected_gene_count:,}"
)
print(
    "Kendall's W: "
    f"{format_number(kendalls_w)}"
)
print(
    "Three-way same direction: "
    f"{unanimous_direction_count:,}/"
    f"{selected_gene_count:,} "
    f"({format_percentage(unanimous_direction_count, selected_gene_count)})"
)