from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np
import numpy.typing as npt
import polars as pl
from scipy.stats import hypergeom, spearmanr

from src.common.paths import (
    resolve_dir_from_root,
    try_find_file_in_project,
)


##########################################################################################
#                                       SETTINGS
##########################################################################################

TAIL_GENE_COUNT: int = 500

analysis_folder: Path = resolve_dir_from_root(
    "data/mtec_low/analysis/affy_vs_rna_seq"
)

analysis_folder.mkdir(
    parents=True,
    exist_ok=True,
)

output_file: Path = (
    analysis_folder
    / f"top_{TAIL_GENE_COUNT}_up_down_agreement.txt"
)

analysed_file: Path = try_find_file_in_project(
    "union_genes.csv"
)

EXPECTED_DATASETS: set[str] = {
    "mouse_gse205410_mtec_low",
    "mouse_gse2585_mtec_low",
    "mouse_gse65617_mtec_low",
}

REQUIRED_COLUMNS: set[str] = {
    "gene_symbol",
    "log2_fold_change",
    "dataset_id",
}


##########################################################################################
#                                     RESULT CLASSES
##########################################################################################


@dataclass(frozen=True, slots=True)
class DatasetTails:
    dataset_id: str

    down_genes: frozenset[str]
    up_genes: frozenset[str]

    negative_gene_count: int
    positive_gene_count: int
    zero_gene_count: int

    down_cutoff: float
    up_cutoff: float

    down_median_effect: float
    up_median_effect: float


@dataclass(frozen=True, slots=True)
class TailReplicationResult:
    direction: str

    same_direction_count: int
    opposite_direction_count: int
    zero_effect_count: int

    same_top_tail_count: int
    opposite_top_tail_count: int

    comparison_median_effect: float


@dataclass(frozen=True, slots=True)
class ReferenceComparisonResult:
    reference_dataset: str
    comparison_dataset: str

    down_result: TailReplicationResult
    up_result: TailReplicationResult

    combined_same_direction_count: int
    combined_opposite_direction_count: int
    combined_zero_effect_count: int

    combined_spearman_rho: float
    combined_spearman_p_value: float


@dataclass(frozen=True, slots=True)
class PairwiseTailOverlapResult:
    first_dataset: str
    second_dataset: str

    down_overlap: int
    down_jaccard: float
    down_enrichment: float
    down_p_value: float

    up_overlap: int
    up_jaccard: float
    up_enrichment: float
    up_p_value: float

    first_down_second_up: int
    first_up_second_down: int


@dataclass(frozen=True, slots=True)
class ThreeWayReferenceResult:
    reference_dataset: str

    down_same_direction_in_all_three: int
    up_same_direction_in_all_three: int

    down_opposite_in_at_least_one: int
    up_opposite_in_at_least_one: int


##########################################################################################
#                                       FORMATTERS
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


##########################################################################################
#                                       DATA HELPERS
##########################################################################################


def get_complete_shared_genes(
    data: pl.DataFrame,
    *,
    dataset_ids: list[str],
) -> pl.DataFrame:
    """
    Retain genes with exactly one finite fold-change value in every dataset.
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


def get_gene_effect_frame(
    wide_frame: pl.DataFrame,
    *,
    genes: frozenset[str] | set[str],
    dataset_ids: list[str],
) -> pl.DataFrame:
    """
    Extract selected genes in deterministic gene-symbol order.
    """
    gene_list: list[str] = sorted(genes)

    selected_frame: pl.DataFrame = (
        wide_frame
        .filter(
            pl.col("gene_symbol").is_in(gene_list)
        )
        .select(
            "gene_symbol",
            *dataset_ids,
        )
        .sort("gene_symbol")
    )

    if selected_frame.height != len(gene_list):
        raise ValueError(
            f"Requested {len(gene_list):,} genes, but extracted "
            f"{selected_frame.height:,}."
        )

    return selected_frame


##########################################################################################
#                                     SELECT TAILS
##########################################################################################


def select_dataset_tails(
    wide_frame: pl.DataFrame,
    *,
    dataset_id: str,
    tail_gene_count: int,
) -> DatasetTails:
    """
    Independently select the most negative and most positive genes.

    Only negative values are eligible for the downregulated list.
    Only positive values are eligible for the upregulated list.
    """
    dataset_frame: pl.DataFrame = (
        wide_frame
        .select(
            "gene_symbol",
            pl.col(dataset_id).alias("effect"),
        )
    )

    negative_frame: pl.DataFrame = dataset_frame.filter(
        pl.col("effect") < 0.0
    )

    positive_frame: pl.DataFrame = dataset_frame.filter(
        pl.col("effect") > 0.0
    )

    zero_frame: pl.DataFrame = dataset_frame.filter(
        pl.col("effect") == 0.0
    )

    if negative_frame.height < tail_gene_count:
        raise ValueError(
            f"Dataset '{dataset_id}' contains only "
            f"{negative_frame.height:,} negative genes, but "
            f"{tail_gene_count:,} were requested."
        )

    if positive_frame.height < tail_gene_count:
        raise ValueError(
            f"Dataset '{dataset_id}' contains only "
            f"{positive_frame.height:,} positive genes, but "
            f"{tail_gene_count:,} were requested."
        )

    down_frame: pl.DataFrame = (
        negative_frame
        .sort(
            by=["effect", "gene_symbol"],
            descending=[False, False],
        )
        .head(tail_gene_count)
    )

    up_frame: pl.DataFrame = (
        positive_frame
        .sort(
            by=["effect", "gene_symbol"],
            descending=[True, False],
        )
        .head(tail_gene_count)
    )

    down_genes: frozenset[str] = frozenset(
        str(value)
        for value in down_frame
        .get_column("gene_symbol")
        .to_list()
    )

    up_genes: frozenset[str] = frozenset(
        str(value)
        for value in up_frame
        .get_column("gene_symbol")
        .to_list()
    )

    if len(down_genes) != tail_gene_count:
        raise ValueError(
            f"Expected {tail_gene_count:,} unique downregulated genes "
            f"for '{dataset_id}', but found {len(down_genes):,}."
        )

    if len(up_genes) != tail_gene_count:
        raise ValueError(
            f"Expected {tail_gene_count:,} unique upregulated genes "
            f"for '{dataset_id}', but found {len(up_genes):,}."
        )

    if down_genes & up_genes:
        raise ValueError(
            f"Dataset '{dataset_id}' has genes in both tails."
        )

    down_values: npt.NDArray[np.float64] = (
        down_frame
        .get_column("effect")
        .to_numpy()
        .astype(np.float64)
    )

    up_values: npt.NDArray[np.float64] = (
        up_frame
        .get_column("effect")
        .to_numpy()
        .astype(np.float64)
    )

    return DatasetTails(
        dataset_id=dataset_id,
        down_genes=down_genes,
        up_genes=up_genes,
        negative_gene_count=negative_frame.height,
        positive_gene_count=positive_frame.height,
        zero_gene_count=zero_frame.height,
        # Closest-to-zero selected negative value.
        down_cutoff=float(np.max(down_values)),
        # Closest-to-zero selected positive value.
        up_cutoff=float(np.min(up_values)),
        down_median_effect=float(np.median(down_values)),
        up_median_effect=float(np.median(up_values)),
    )


##########################################################################################
#                                 TAIL REPLICATION
##########################################################################################


def evaluate_tail_replication(
    wide_frame: pl.DataFrame,
    *,
    reference_dataset: str,
    comparison_dataset: str,
    reference_genes: frozenset[str],
    reference_direction: str,
    comparison_tails: DatasetTails,
) -> TailReplicationResult:
    """
    Evaluate one reference tail in another dataset.

    For example, take the 500 most decreased genes from GSE205410 and
    determine how many are also negative, top-down, positive, or top-up
    in GSE2585.
    """
    if reference_direction not in {"down", "up"}:
        raise ValueError(
            "reference_direction must be 'down' or 'up'."
        )

    selected_frame: pl.DataFrame = get_gene_effect_frame(
        wide_frame,
        genes=reference_genes,
        dataset_ids=[
            reference_dataset,
            comparison_dataset,
        ],
    )

    comparison_values: npt.NDArray[np.float64] = (
        selected_frame
        .get_column(comparison_dataset)
        .to_numpy()
        .astype(np.float64)
    )

    if reference_direction == "down":
        same_direction_count: int = int(
            np.sum(comparison_values < 0.0)
        )

        opposite_direction_count: int = int(
            np.sum(comparison_values > 0.0)
        )

        same_top_tail_count: int = len(
            reference_genes
            & comparison_tails.down_genes
        )

        opposite_top_tail_count: int = len(
            reference_genes
            & comparison_tails.up_genes
        )

    else:
        same_direction_count = int(
            np.sum(comparison_values > 0.0)
        )

        opposite_direction_count = int(
            np.sum(comparison_values < 0.0)
        )

        same_top_tail_count = len(
            reference_genes
            & comparison_tails.up_genes
        )

        opposite_top_tail_count = len(
            reference_genes
            & comparison_tails.down_genes
        )

    zero_effect_count: int = int(
        np.sum(comparison_values == 0.0)
    )

    return TailReplicationResult(
        direction=reference_direction,
        same_direction_count=same_direction_count,
        opposite_direction_count=opposite_direction_count,
        zero_effect_count=zero_effect_count,
        same_top_tail_count=same_top_tail_count,
        opposite_top_tail_count=opposite_top_tail_count,
        comparison_median_effect=float(
            np.median(comparison_values)
        ),
    )


def evaluate_reference_comparison(
    wide_frame: pl.DataFrame,
    *,
    reference_tails: DatasetTails,
    comparison_tails: DatasetTails,
) -> ReferenceComparisonResult:
    """
    Evaluate one dataset's 500 down + 500 up genes in another dataset.
    """
    reference_dataset: str = (
        reference_tails.dataset_id
    )

    comparison_dataset: str = (
        comparison_tails.dataset_id
    )

    down_result: TailReplicationResult = (
        evaluate_tail_replication(
            wide_frame,
            reference_dataset=reference_dataset,
            comparison_dataset=comparison_dataset,
            reference_genes=reference_tails.down_genes,
            reference_direction="down",
            comparison_tails=comparison_tails,
        )
    )

    up_result: TailReplicationResult = (
        evaluate_tail_replication(
            wide_frame,
            reference_dataset=reference_dataset,
            comparison_dataset=comparison_dataset,
            reference_genes=reference_tails.up_genes,
            reference_direction="up",
            comparison_tails=comparison_tails,
        )
    )

    combined_genes: frozenset[str] = (
        reference_tails.down_genes
        | reference_tails.up_genes
    )

    combined_frame: pl.DataFrame = get_gene_effect_frame(
        wide_frame,
        genes=combined_genes,
        dataset_ids=[
            reference_dataset,
            comparison_dataset,
        ],
    )

    reference_values: npt.NDArray[np.float64] = (
        combined_frame
        .get_column(reference_dataset)
        .to_numpy()
        .astype(np.float64)
    )

    comparison_values: npt.NDArray[np.float64] = (
        combined_frame
        .get_column(comparison_dataset)
        .to_numpy()
        .astype(np.float64)
    )

    reference_signs: npt.NDArray[np.float64] = np.sign(
        reference_values
    )

    comparison_signs: npt.NDArray[np.float64] = np.sign(
        comparison_values
    )

    combined_same_direction_count: int = int(
        np.sum(
            (reference_signs == comparison_signs)
            & (comparison_signs != 0.0)
        )
    )

    combined_opposite_direction_count: int = int(
        np.sum(
            reference_signs
            == -comparison_signs
        )
    )

    combined_zero_effect_count: int = int(
        np.sum(
            comparison_signs == 0.0
        )
    )

    spearman_result = spearmanr(
        reference_values,
        comparison_values,
        nan_policy="raise",
    )

    return ReferenceComparisonResult(
        reference_dataset=reference_dataset,
        comparison_dataset=comparison_dataset,
        down_result=down_result,
        up_result=up_result,
        combined_same_direction_count=(
            combined_same_direction_count
        ),
        combined_opposite_direction_count=(
            combined_opposite_direction_count
        ),
        combined_zero_effect_count=(
            combined_zero_effect_count
        ),
        combined_spearman_rho=float(
            spearman_result.statistic
        ),
        combined_spearman_p_value=float(
            spearman_result.pvalue
        ),
    )


##########################################################################################
#                                  PAIRWISE OVERLAPS
##########################################################################################


def calculate_jaccard(
    first_set: frozenset[str],
    second_set: frozenset[str],
) -> float:
    union_size: int = len(
        first_set | second_set
    )

    if union_size == 0:
        return float("nan")

    return (
        len(first_set & second_set)
        / union_size
    )


def calculate_overlap_p_value(
    *,
    overlap: int,
    universe_size: int,
    first_set_size: int,
    second_set_size: int,
) -> float:
    """
    Hypergeometric upper-tail probability.

    This asks how likely an overlap at least this large would be if both
    gene lists were randomly selected from the shared-gene universe.
    """
    return float(
        hypergeom.sf(
            overlap - 1,
            universe_size,
            first_set_size,
            second_set_size,
        )
    )


def calculate_pairwise_tail_overlap(
    first_tails: DatasetTails,
    second_tails: DatasetTails,
    *,
    universe_size: int,
) -> PairwiseTailOverlapResult:
    down_overlap: int = len(
        first_tails.down_genes
        & second_tails.down_genes
    )

    up_overlap: int = len(
        first_tails.up_genes
        & second_tails.up_genes
    )

    first_down_second_up: int = len(
        first_tails.down_genes
        & second_tails.up_genes
    )

    first_up_second_down: int = len(
        first_tails.up_genes
        & second_tails.down_genes
    )

    expected_overlap: float = (
        len(first_tails.down_genes)
        * len(second_tails.down_genes)
        / universe_size
    )

    down_enrichment: float = (
        down_overlap / expected_overlap
        if expected_overlap > 0
        else float("nan")
    )

    up_enrichment: float = (
        up_overlap / expected_overlap
        if expected_overlap > 0
        else float("nan")
    )

    return PairwiseTailOverlapResult(
        first_dataset=first_tails.dataset_id,
        second_dataset=second_tails.dataset_id,
        down_overlap=down_overlap,
        down_jaccard=calculate_jaccard(
            first_tails.down_genes,
            second_tails.down_genes,
        ),
        down_enrichment=down_enrichment,
        down_p_value=calculate_overlap_p_value(
            overlap=down_overlap,
            universe_size=universe_size,
            first_set_size=len(
                first_tails.down_genes
            ),
            second_set_size=len(
                second_tails.down_genes
            ),
        ),
        up_overlap=up_overlap,
        up_jaccard=calculate_jaccard(
            first_tails.up_genes,
            second_tails.up_genes,
        ),
        up_enrichment=up_enrichment,
        up_p_value=calculate_overlap_p_value(
            overlap=up_overlap,
            universe_size=universe_size,
            first_set_size=len(
                first_tails.up_genes
            ),
            second_set_size=len(
                second_tails.up_genes
            ),
        ),
        first_down_second_up=(
            first_down_second_up
        ),
        first_up_second_down=(
            first_up_second_down
        ),
    )


##########################################################################################
#                                THREE-WAY REPLICATION
##########################################################################################


def evaluate_three_way_reference(
    wide_frame: pl.DataFrame,
    *,
    reference_tails: DatasetTails,
    dataset_ids: list[str],
) -> ThreeWayReferenceResult:
    other_datasets: list[str] = [
        dataset_id
        for dataset_id in dataset_ids
        if dataset_id != reference_tails.dataset_id
    ]

    if len(other_datasets) != 2:
        raise ValueError(
            "Three-way reference evaluation requires exactly "
            "three datasets."
        )

    first_other: str = other_datasets[0]
    second_other: str = other_datasets[1]

    down_frame: pl.DataFrame = get_gene_effect_frame(
        wide_frame,
        genes=reference_tails.down_genes,
        dataset_ids=[
            first_other,
            second_other,
        ],
    )

    down_matrix: npt.NDArray[np.float64] = (
        down_frame
        .select(
            first_other,
            second_other,
        )
        .to_numpy()
        .astype(np.float64)
    )

    down_same_direction_in_all_three: int = int(
        np.sum(
            np.all(
                down_matrix < 0.0,
                axis=1,
            )
        )
    )

    down_opposite_in_at_least_one: int = int(
        np.sum(
            np.any(
                down_matrix > 0.0,
                axis=1,
            )
        )
    )

    up_frame: pl.DataFrame = get_gene_effect_frame(
        wide_frame,
        genes=reference_tails.up_genes,
        dataset_ids=[
            first_other,
            second_other,
        ],
    )

    up_matrix: npt.NDArray[np.float64] = (
        up_frame
        .select(
            first_other,
            second_other,
        )
        .to_numpy()
        .astype(np.float64)
    )

    up_same_direction_in_all_three: int = int(
        np.sum(
            np.all(
                up_matrix > 0.0,
                axis=1,
            )
        )
    )

    up_opposite_in_at_least_one: int = int(
        np.sum(
            np.any(
                up_matrix < 0.0,
                axis=1,
            )
        )
    )

    return ThreeWayReferenceResult(
        reference_dataset=(
            reference_tails.dataset_id
        ),
        down_same_direction_in_all_three=(
            down_same_direction_in_all_three
        ),
        up_same_direction_in_all_three=(
            up_same_direction_in_all_three
        ),
        down_opposite_in_at_least_one=(
            down_opposite_in_at_least_one
        ),
        up_opposite_in_at_least_one=(
            up_opposite_in_at_least_one
        ),
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
    REQUIRED_COLUMNS
    - set(source_df.columns)
)

if missing_columns:
    raise ValueError(
        "Input file is missing required columns: "
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
#                                       VALIDATE
##########################################################################################

observed_datasets: set[str] = set(
    source_df
    .get_column("dataset_id")
    .unique()
    .to_list()
)

if observed_datasets != EXPECTED_DATASETS:
    missing_datasets: set[str] = (
        EXPECTED_DATASETS
        - observed_datasets
    )

    unexpected_datasets: set[str] = (
        observed_datasets
        - EXPECTED_DATASETS
    )

    raise ValueError(
        "Dataset IDs do not match expected datasets.\n"
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
    str(row["dataset_id"]): int(
        row["n_genes"]
    )
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
#                                 COMPLETE SHARED GENES
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

effect_matrix: npt.NDArray[np.float64] = (
    wide_df
    .select(dataset_ids)
    .to_numpy()
    .astype(np.float64)
)

if not np.isfinite(effect_matrix).all():
    raise ValueError(
        "The complete shared-gene matrix contains "
        "non-finite values."
    )


##########################################################################################
#                                  SELECT SIX LISTS
##########################################################################################

tails_by_dataset: dict[str, DatasetTails] = {
    dataset_id: select_dataset_tails(
        wide_df,
        dataset_id=dataset_id,
        tail_gene_count=TAIL_GENE_COUNT,
    )
    for dataset_id in dataset_ids
}


##########################################################################################
#                                PAIRWISE LIST OVERLAP
##########################################################################################

pairwise_overlap_results: list[
    PairwiseTailOverlapResult
] = []

for first_dataset, second_dataset in combinations(
    dataset_ids,
    2,
):
    pairwise_overlap_results.append(
        calculate_pairwise_tail_overlap(
            tails_by_dataset[first_dataset],
            tails_by_dataset[second_dataset],
            universe_size=shared_gene_count,
        )
    )


##########################################################################################
#                           REFERENCE-SPECIFIC REPLICATION
##########################################################################################

reference_comparison_results: list[
    ReferenceComparisonResult
] = []

for reference_dataset in dataset_ids:
    for comparison_dataset in dataset_ids:
        if reference_dataset == comparison_dataset:
            continue

        reference_comparison_results.append(
            evaluate_reference_comparison(
                wide_df,
                reference_tails=(
                    tails_by_dataset[
                        reference_dataset
                    ]
                ),
                comparison_tails=(
                    tails_by_dataset[
                        comparison_dataset
                    ]
                ),
            )
        )


##########################################################################################
#                               THREE-WAY AGREEMENT
##########################################################################################

common_top_down_genes: set[str] = set(
    tails_by_dataset[
        dataset_ids[0]
    ].down_genes
)

common_top_up_genes: set[str] = set(
    tails_by_dataset[
        dataset_ids[0]
    ].up_genes
)

for dataset_id in dataset_ids[1:]:
    common_top_down_genes.intersection_update(
        tails_by_dataset[
            dataset_id
        ].down_genes
    )

    common_top_up_genes.intersection_update(
        tails_by_dataset[
            dataset_id
        ].up_genes
    )

three_way_reference_results: list[
    ThreeWayReferenceResult
] = [
    evaluate_three_way_reference(
        wide_df,
        reference_tails=(
            tails_by_dataset[dataset_id]
        ),
        dataset_ids=dataset_ids,
    )
    for dataset_id in dataset_ids
]

expected_pairwise_overlap: float = (
    TAIL_GENE_COUNT**2
    / shared_gene_count
)

expected_three_way_overlap: float = (
    TAIL_GENE_COUNT**3
    / shared_gene_count**2
)


##########################################################################################
#                                       BUILD REPORT
##########################################################################################

report_lines: list[str] = [
    "MTEC-LOW TOP-500 UP/DOWN AGREEMENT REPORT",
    "=" * 80,
    "",
    "Analysis objective",
    "-" * 80,
    (
        f"Each dataset independently contributes its "
        f"{TAIL_GENE_COUNT:,} most decreased genes and "
        f"{TAIL_GENE_COUNT:,} most increased genes."
    ),
    (
        "The other datasets are then evaluated against each reference "
        "list. No consensus gene selection is performed."
    ),
    (
        "This allows effects detected by only one dataset to remain "
        "visible rather than requiring all datasets to show a strong effect."
    ),
    "",
    "Direction definition",
    "-" * 80,
    (
        "Decreased: negative KO-versus-control log2 fold change."
    ),
    (
        "Increased: positive KO-versus-control log2 fold change."
    ),
    (
        "The analysis uses fold-change direction and within-dataset "
        "ranking. Raw fold-change magnitudes are not compared between "
        "platforms."
    ),
    "",
    "Input",
    "-" * 80,
    f"Source file: {analysed_file}",
    f"Input rows: {source_df.height:,}",
    f"Unique genes in input: {total_unique_genes:,}",
    (
        "Genes with one finite effect in every dataset: "
        f"{shared_gene_count:,}"
    ),
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
        "Independent tail selection",
        "-" * 80,
    ]
)

for dataset_id in dataset_ids:
    tails: DatasetTails = (
        tails_by_dataset[dataset_id]
    )

    report_lines.extend(
        [
            dataset_id,
            (
                "  Available negative genes: "
                f"{tails.negative_gene_count:,}"
            ),
            (
                "  Available positive genes: "
                f"{tails.positive_gene_count:,}"
            ),
            (
                "  Exact-zero genes: "
                f"{tails.zero_gene_count:,}"
            ),
            (
                f"  Selected decreased genes: "
                f"{len(tails.down_genes):,}"
            ),
            (
                "  Decreased-list boundary log2FC: "
                f"{format_number(tails.down_cutoff)}"
            ),
            (
                "  Median selected decreased log2FC: "
                f"{format_number(tails.down_median_effect)}"
            ),
            (
                f"  Selected increased genes: "
                f"{len(tails.up_genes):,}"
            ),
            (
                "  Increased-list boundary log2FC: "
                f"{format_number(tails.up_cutoff)}"
            ),
            (
                "  Median selected increased log2FC: "
                f"{format_number(tails.up_median_effect)}"
            ),
            "",
        ]
    )

report_lines.extend(
    [
        "Pairwise overlap of independently selected top lists",
        "-" * 80,
        (
            "The overlap count asks how many genes appear in the same "
            "top-500 directional list in both datasets."
        ),
        (
            f"With {shared_gene_count:,} shared genes, two random "
            f"500-gene lists would overlap by approximately "
            f"{expected_pairwise_overlap:.1f} genes."
        ),
        "",
    ]
)

for result in pairwise_overlap_results:
    same_tail_total: int = (
        result.down_overlap
        + result.up_overlap
    )

    opposite_tail_total: int = (
        result.first_down_second_up
        + result.first_up_second_down
    )

    report_lines.extend(
        [
            result.first_dataset,
            f"versus {result.second_dataset}",
            (
                "  Shared top-decreased genes: "
                f"{result.down_overlap:,}/"
                f"{TAIL_GENE_COUNT:,} "
                f"({format_percentage(result.down_overlap, TAIL_GENE_COUNT)})"
            ),
            (
                "    Jaccard similarity: "
                f"{format_number(result.down_jaccard)}"
            ),
            (
                "    Enrichment over random overlap: "
                f"{format_number(result.down_enrichment)}-fold"
            ),
            (
                "    Hypergeometric p-value: "
                f"{format_p_value(result.down_p_value)}"
            ),
            (
                "  Shared top-increased genes: "
                f"{result.up_overlap:,}/"
                f"{TAIL_GENE_COUNT:,} "
                f"({format_percentage(result.up_overlap, TAIL_GENE_COUNT)})"
            ),
            (
                "    Jaccard similarity: "
                f"{format_number(result.up_jaccard)}"
            ),
            (
                "    Enrichment over random overlap: "
                f"{format_number(result.up_enrichment)}-fold"
            ),
            (
                "    Hypergeometric p-value: "
                f"{format_p_value(result.up_p_value)}"
            ),
            (
                "  Same strong tail in both datasets: "
                f"{same_tail_total:,}/"
                f"{2 * TAIL_GENE_COUNT:,} "
                f"({format_percentage(same_tail_total, 2 * TAIL_GENE_COUNT)})"
            ),
            (
                "  First dataset top-decreased but second top-increased: "
                f"{result.first_down_second_up:,}"
            ),
            (
                "  First dataset top-increased but second top-decreased: "
                f"{result.first_up_second_down:,}"
            ),
            (
                "  Strong opposite-tail contradictions: "
                f"{opposite_tail_total:,}/"
                f"{2 * TAIL_GENE_COUNT:,} "
                f"({format_percentage(opposite_tail_total, 2 * TAIL_GENE_COUNT)})"
            ),
            "",
        ]
    )

report_lines.extend(
    [
        "Reference-specific replication",
        "-" * 80,
        (
            "Each section selects 500 decreased and 500 increased genes "
            "from the named reference dataset only."
        ),
        (
            "The comparison dataset is not involved in selecting those genes."
        ),
        "",
    ]
)

for result in reference_comparison_results:
    combined_gene_count: int = (
        2 * TAIL_GENE_COUNT
    )

    report_lines.extend(
        [
            f"Reference: {result.reference_dataset}",
            f"Evaluated in: {result.comparison_dataset}",
            "",
            (
                f"  Reference top-{TAIL_GENE_COUNT:,} decreased genes"
            ),
            (
                "    Negative in comparison dataset: "
                f"{result.down_result.same_direction_count:,}/"
                f"{TAIL_GENE_COUNT:,} "
                f"({format_percentage(result.down_result.same_direction_count, TAIL_GENE_COUNT)})"
            ),
            (
                "    Positive in comparison dataset: "
                f"{result.down_result.opposite_direction_count:,}/"
                f"{TAIL_GENE_COUNT:,} "
                f"({format_percentage(result.down_result.opposite_direction_count, TAIL_GENE_COUNT)})"
            ),
            (
                "    Also comparison top-decreased: "
                f"{result.down_result.same_top_tail_count:,}/"
                f"{TAIL_GENE_COUNT:,} "
                f"({format_percentage(result.down_result.same_top_tail_count, TAIL_GENE_COUNT)})"
            ),
            (
                "    Comparison top-increased contradiction: "
                f"{result.down_result.opposite_top_tail_count:,}/"
                f"{TAIL_GENE_COUNT:,} "
                f"({format_percentage(result.down_result.opposite_top_tail_count, TAIL_GENE_COUNT)})"
            ),
            (
                "    Median comparison log2FC: "
                f"{format_number(result.down_result.comparison_median_effect)}"
            ),
            "",
            (
                f"  Reference top-{TAIL_GENE_COUNT:,} increased genes"
            ),
            (
                "    Positive in comparison dataset: "
                f"{result.up_result.same_direction_count:,}/"
                f"{TAIL_GENE_COUNT:,} "
                f"({format_percentage(result.up_result.same_direction_count, TAIL_GENE_COUNT)})"
            ),
            (
                "    Negative in comparison dataset: "
                f"{result.up_result.opposite_direction_count:,}/"
                f"{TAIL_GENE_COUNT:,} "
                f"({format_percentage(result.up_result.opposite_direction_count, TAIL_GENE_COUNT)})"
            ),
            (
                "    Also comparison top-increased: "
                f"{result.up_result.same_top_tail_count:,}/"
                f"{TAIL_GENE_COUNT:,} "
                f"({format_percentage(result.up_result.same_top_tail_count, TAIL_GENE_COUNT)})"
            ),
            (
                "    Comparison top-decreased contradiction: "
                f"{result.up_result.opposite_top_tail_count:,}/"
                f"{TAIL_GENE_COUNT:,} "
                f"({format_percentage(result.up_result.opposite_top_tail_count, TAIL_GENE_COUNT)})"
            ),
            (
                "    Median comparison log2FC: "
                f"{format_number(result.up_result.comparison_median_effect)}"
            ),
            "",
            (
                f"  Combined reference list: "
                f"{combined_gene_count:,} genes"
            ),
            (
                "    Same direction in comparison dataset: "
                f"{result.combined_same_direction_count:,}/"
                f"{combined_gene_count:,} "
                f"({format_percentage(result.combined_same_direction_count, combined_gene_count)})"
            ),
            (
                "    Opposite direction in comparison dataset: "
                f"{result.combined_opposite_direction_count:,}/"
                f"{combined_gene_count:,} "
                f"({format_percentage(result.combined_opposite_direction_count, combined_gene_count)})"
            ),
            (
                "    Exact zero in comparison dataset: "
                f"{result.combined_zero_effect_count:,}/"
                f"{combined_gene_count:,} "
                f"({format_percentage(result.combined_zero_effect_count, combined_gene_count)})"
            ),
            (
                "    Spearman rho across the combined reference list: "
                f"{format_number(result.combined_spearman_rho)}"
            ),
            (
                "    Descriptive Spearman p-value: "
                f"{format_p_value(result.combined_spearman_p_value)}"
            ),
            "",
        ]
    )

report_lines.extend(
    [
        "Three-dataset agreement",
        "-" * 80,
        (
            "Common top-decreased genes in all three datasets: "
            f"{len(common_top_down_genes):,}/"
            f"{TAIL_GENE_COUNT:,} "
            f"({format_percentage(len(common_top_down_genes), TAIL_GENE_COUNT)})"
        ),
        (
            "Common top-increased genes in all three datasets: "
            f"{len(common_top_up_genes):,}/"
            f"{TAIL_GENE_COUNT:,} "
            f"({format_percentage(len(common_top_up_genes), TAIL_GENE_COUNT)})"
        ),
        (
            "Expected random overlap of three independent 500-gene "
            f"lists: approximately {expected_three_way_overlap:.2f} genes."
        ),
        "",
        (
            "The following values ask whether each reference dataset's "
            "top genes at least retain the same sign in both other datasets."
        ),
        "",
    ]
)

for result in three_way_reference_results:
    report_lines.extend(
        [
            f"Reference: {result.reference_dataset}",
            (
                "  Reference decreased genes negative in both "
                "other datasets: "
                f"{result.down_same_direction_in_all_three:,}/"
                f"{TAIL_GENE_COUNT:,} "
                f"({format_percentage(result.down_same_direction_in_all_three, TAIL_GENE_COUNT)})"
            ),
            (
                "  Reference decreased genes positive in at least "
                "one other dataset: "
                f"{result.down_opposite_in_at_least_one:,}/"
                f"{TAIL_GENE_COUNT:,} "
                f"({format_percentage(result.down_opposite_in_at_least_one, TAIL_GENE_COUNT)})"
            ),
            (
                "  Reference increased genes positive in both "
                "other datasets: "
                f"{result.up_same_direction_in_all_three:,}/"
                f"{TAIL_GENE_COUNT:,} "
                f"({format_percentage(result.up_same_direction_in_all_three, TAIL_GENE_COUNT)})"
            ),
            (
                "  Reference increased genes negative in at least "
                "one other dataset: "
                f"{result.up_opposite_in_at_least_one:,}/"
                f"{TAIL_GENE_COUNT:,} "
                f"({format_percentage(result.up_opposite_in_at_least_one, TAIL_GENE_COUNT)})"
            ),
            "",
        ]
    )

report_lines.extend(
    [
        "Interpretation notes",
        "-" * 80,
        (
            "Same sign is the broadest replication criterion. For example, "
            "a reference-decreased gene counts as replicated when its "
            "comparison fold change is any negative value."
        ),
        (
            "Same top tail is stricter. It requires the gene to be among "
            "the 500 strongest effects in the same direction in both datasets."
        ),
        (
            "Opposite top tail is a strong contradiction because the gene "
            "is among the strongest decreased genes in one dataset and the "
            "strongest increased genes in another."
        ),
        (
            "Spearman rho evaluates whether the combined 500 decreased and "
            "500 increased reference genes have similar ordering in the "
            "comparison dataset. It does not require equal fold-change scales."
        ),
        (
            "The hypergeometric p-value tests whether the top-list overlap "
            "is larger than expected from random 500-gene lists drawn from "
            "the shared-gene universe."
        ),
        (
            "Because genes are selected using observed fold changes, all "
            "reported p-values should be treated as descriptive rather than "
            "independent confirmatory tests."
        ),
        (
            "Adjusted differential-expression p-values are not used for "
            "selection because sample sizes and statistical power differ "
            "between datasets."
        ),
        (
            "A small fold change means little KO-versus-control difference; "
            "it does not necessarily mean low absolute gene expression."
        ),
        (
            "Genes tied exactly at a top-list boundary are resolved "
            "deterministically by gene symbol so every list contains "
            f"exactly {TAIL_GENE_COUNT:,} genes."
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
print(f"Shared-gene universe: {shared_gene_count:,}")
print(
    f"Selected per dataset: "
    f"{TAIL_GENE_COUNT:,} decreased + "
    f"{TAIL_GENE_COUNT:,} increased"
)
print(
    "Common top-decreased genes in all three: "
    f"{len(common_top_down_genes):,}"
)
print(
    "Common top-increased genes in all three: "
    f"{len(common_top_up_genes):,}"
)