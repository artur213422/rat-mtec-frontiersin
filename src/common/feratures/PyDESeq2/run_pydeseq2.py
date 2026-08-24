from __future__ import annotations

import numpy as np
import pandas as pd
import polars as pl
from numpy.typing import NDArray
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats


def _validate_sample_columns(counts: pl.DataFrame, control_columns: list[str], ko_columns: list[str],) -> list[str]:
    if not control_columns:
        raise ValueError("control_columns cannot be empty")

    if not ko_columns:
        raise ValueError("ko_columns cannot be empty")

    duplicated_columns: set[str] = (
        set(control_columns) & set(ko_columns)
    )

    if duplicated_columns:
        raise ValueError("The following samples occur in both groups: {sorted(duplicated_columns)}")

    sample_columns: list[str] = [
        *control_columns,
        *ko_columns,
    ]

    missing_columns: set[str] = (
        set(sample_columns) - set(counts.columns)
    )

    if missing_columns:
        raise ValueError(f"Missing sample columns: {sorted(missing_columns)}. Available columns: {counts.columns}")

    return sample_columns


def _extract_integer_count_matrix(counts: pl.DataFrame, sample_columns: list[str]) -> NDArray[np.int64]:
    if counts.height == 0:
        raise ValueError("The count matrix contains no rows")

    selected_counts: pl.DataFrame = counts.select(sample_columns)

    try:
        numeric_matrix: NDArray[np.float64] = np.asarray(
            selected_counts.to_numpy(),
            dtype=np.float64,
        )
    except (TypeError, ValueError) as error:
        raise ValueError("All selected sample columns must contain numeric values") from error

    if not np.isfinite(numeric_matrix).all():
        raise ValueError("The count matrix contains null, NaN, or infinite values")

    if (numeric_matrix < 0).any():
        raise ValueError("PyDESeq2 raw counts cannot contain negative values")

    if not np.equal(numeric_matrix, np.floor(numeric_matrix)).all():
        raise ValueError("PyDESeq2 requires raw integer counts. Non-integer values were found.")

    if (numeric_matrix > np.iinfo(np.int64).max).any():
        raise ValueError("Some count values are too large for Int64")

    return numeric_matrix.astype(np.int64, copy=False,)


def _create_pydeseq2_inputs(count_matrix: NDArray[np.int64], control_columns: list[str], ko_columns: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_columns: list[str] = [
        *control_columns,
        *ko_columns,
    ]

    feature_names: list[str] = [
        f"feature_{row_index}"
        for row_index in range(count_matrix.shape[0])
    ]

    # Input is currently genes × samples.
    # PyDESeq2 requires samples × genes.
    counts_pd: pd.DataFrame = pd.DataFrame(
        data=count_matrix.transpose(),
        index=pd.Index(
            sample_columns,
            name="sample",
        ),
        columns=feature_names,
        dtype=np.int64,
    )

    conditions: list[str] = [
        *(["control"] * len(control_columns)),
        *(["ko"] * len(ko_columns)),
    ]

    metadata_pd: pd.DataFrame = pd.DataFrame(
        {
            "condition": pd.Categorical(
                conditions,
                categories=["control", "ko"],
            )
        },
        index=pd.Index(
            sample_columns,
            name="sample",
        ),
    )

    return counts_pd, metadata_pd


def _convert_results_to_polars(results: pd.DataFrame) -> pl.DataFrame:
    results_pd: pd.DataFrame = results.copy()

    # Recover input row order from temporary feature names.
    results_pd["_row_index"] = [
        int(str(feature_name).removeprefix("feature_"))
        for feature_name in results_pd.index
    ]

    results_pd = (
        results_pd
        .sort_values("_row_index")
        .drop(columns="_row_index")
        .reset_index(drop=True)
    )

    result: pl.DataFrame = pl.from_pandas(
        results_pd,
        include_index=False,
    )

    numeric_columns: list[str] = [
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj",
    ]

    return result.with_columns(
        [
            pl.col(column)
            .cast(pl.Float64)
            .fill_nan(None)
            for column in numeric_columns
            if column in result.columns
        ]
    )

def run_pydeseq2(
    counts: pl.DataFrame,
    control_columns: list[str],
    ko_columns: list[str],
    *,
    metadata_columns: list[str] | None = None,
    alpha: float = 0.05,
    n_cpus: int = 1,
) -> pl.DataFrame:
    """
    Run KO-versus-control differential-expression analysis.

    Positive log2FoldChange means higher expression in KO.

    Metadata columns are preserved in the returned DataFrame but are
    not used by PyDESeq2.
    """
    if metadata_columns is None:
        metadata_columns = [
            "gene_symbol",
            "has_clean_ortholog",
            "ortholog_id",
        ]

    if not 0.0 < alpha < 1.0:
        raise ValueError(
            f"alpha must be between 0 and 1, received {alpha}"
        )

    if n_cpus < 1:
        raise ValueError(
            f"n_cpus must be at least 1, received {n_cpus}"
        )

    missing_metadata_columns: set[str] = (
        set(metadata_columns) - set(counts.columns)
    )

    if missing_metadata_columns:
        raise ValueError(
            "Missing metadata columns: "
            f"{sorted(missing_metadata_columns)}. "
            f"Available columns: {counts.columns}"
        )

    sample_columns: list[str] = _validate_sample_columns(
        counts,
        control_columns,
        ko_columns,
    )

    count_matrix: NDArray[np.int64] = (
        _extract_integer_count_matrix(
            counts,
            sample_columns,
        )
    )

    counts_pd, metadata_pd = _create_pydeseq2_inputs(
        count_matrix,
        control_columns,
        ko_columns,
    )

    inference = DefaultInference(
        n_cpus=n_cpus,
    )

    dds = DeseqDataSet(
        counts=counts_pd,
        metadata=metadata_pd,
        design="~condition",
        refit_cooks=True,
        inference=inference,
        quiet=True,
    )

    dds.deseq2()

    statistics = DeseqStats(
        dds,
        contrast=[
            "condition",
            "ko",
            "control",
        ],
        alpha=alpha,
        inference=inference,
        quiet=True,
    )

    statistics.summary()

    result: pl.DataFrame = _convert_results_to_polars(
        statistics.results_df
    )

    if result.height != counts.height:
        raise RuntimeError(
            "PyDESeq2 result row count does not match input row count: "
            f"{result.height} != {counts.height}"
        )

    return pl.concat(
        [
            counts.select(metadata_columns),
            result,
        ],
        how="horizontal",
    )

def prepare_pydeseq2(
    counts: pl.DataFrame,
    control_columns: list[str],
    ko_columns: list[str],
) -> pl.DataFrame:
    """
    Temporarily prepare an RNA-seq matrix for PyDESeq2 by rounding
    fractional count values and casting sample columns to Int64.

    Use only when fractional values were introduced by duplicate-gene
    averaging. Do not use this on RPKM, TPM, FPKM, CPM, or other
    normalized expression matrices.
    """
    sample_columns: list[str] = [
        *control_columns,
        *ko_columns,
    ]

    if not control_columns:
        raise ValueError("control_columns cannot be empty")

    if not ko_columns:
        raise ValueError("ko_columns cannot be empty")

    duplicated_columns: set[str] = (
        set(control_columns) & set(ko_columns)
    )

    if duplicated_columns:
        raise ValueError(
            "Columns appear in both control and KO groups: "
            f"{sorted(duplicated_columns)}"
        )

    missing_columns: set[str] = (
        set(sample_columns) - set(counts.columns)
    )

    if missing_columns:
        raise ValueError(
            f"Missing sample columns: {sorted(missing_columns)}"
        )

    null_columns: list[str] = [
        column
        for column in sample_columns
        if counts.get_column(column).null_count() > 0
    ]

    if null_columns:
        raise ValueError(
            f"Null values found in sample columns: {null_columns}"
        )

    has_negative_counts: bool = counts.select(
        pl.any_horizontal(
            [
                pl.col(column) < 0
                for column in sample_columns
            ]
        ).any()
    ).item()

    if has_negative_counts:
        raise ValueError("Negative count values were found")

    return counts.with_columns(
        pl.col(sample_columns)
        .round(0)
        .cast(pl.Int64)
    )