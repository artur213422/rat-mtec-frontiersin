from typing import Literal

import numpy as np
import polars as pl


SelectionRule = Literal[
    "largest_j",
    "smallest_j",
]

ExpressionScale = Literal[
    "log2",
    "linear",
]


def perform_dynamic_step(
    atlas_df: pl.DataFrame,
    *,
    probe_id_column: str = "ID_REF",
    metadata_columns: list[str] | None = None,
    maximum_selected_groups: int = 5,
    minimum_expression: float = 50.0,
    selection_rule: SelectionRule = "largest_j",
    expression_scale: ExpressionScale = "log2",
) -> pl.DataFrame:
    """
    Apply the Sansom dynamic-step classifier to a probe-by-group matrix.

    The function is pure:
    - it does not read or write files;
    - it does not mutate atlas_df;
    - it does not print;
    - identical inputs produce identical outputs.

    Parameters
    ----------
    atlas_df
        One row per probe. Contains ID_REF and physiological-group
        expression columns.

    probe_id_column
        Name of the probe identifier column.

    maximum_selected_groups
        Maximum number of high-expression groups considered restricted.

    minimum_expression
        Minimum linear-scale expression required at the selected boundary.

    selection_rule
        Rule used when multiple boundaries qualify.

    expression_scale
        Scale of the expression columns. Dynamic-step calculations are
        performed on the linear scale.

    Returns
    -------
    pl.DataFrame
        One row per probe containing its tissue-restriction call.
    """

    if metadata_columns is None:
        metadata_columns = []

    duplicated_metadata_columns: list[str] = sorted(
        {
            column
            for column in metadata_columns
            if metadata_columns.count(column) > 1
        }
    )

    if duplicated_metadata_columns:
        raise ValueError(f"Duplicated metadata columns: {duplicated_metadata_columns}")

    missing_metadata_columns: list[str] = [
        column
        for column in metadata_columns
        if column not in atlas_df.columns
    ]

    if missing_metadata_columns:
        raise ValueError(f"Metadata columns are missing from the atlas: {missing_metadata_columns}")

    if probe_id_column in metadata_columns:
        raise ValueError(f"{probe_id_column} should not also be listed in metadata_columns")

    if probe_id_column not in atlas_df.columns:
        raise ValueError(f"Required probe ID column is missing: {probe_id_column}")

    if maximum_selected_groups < 1:
        raise ValueError("maximum_selected_groups must be at least 1")

    if minimum_expression < 0:
        raise ValueError("minimum_expression cannot be negative")

    if selection_rule not in {
        "largest_j",
        "smallest_j",
    }:
        raise ValueError(f"Unsupported selection rule: {selection_rule}")

    if expression_scale not in {
        "log2",
        "linear",
    }:
        raise ValueError(f"Unsupported expression scale: {expression_scale}")

    excluded_columns: set[str] = {
        probe_id_column,
        *metadata_columns,
    }

    group_columns: list[str] = [
        column
        for column in atlas_df.columns
        if column not in excluded_columns
    ]

    if len(group_columns) <= maximum_selected_groups:
        raise ValueError(f"Found {len(group_columns)} physiological groups, but at least {maximum_selected_groups + 1} are ""required")

    duplicated_probe_count: int = (
        atlas_df
        .filter(pl.col(probe_id_column).is_duplicated())
        .height
    )

    if duplicated_probe_count > 0:
        raise ValueError(f"Found {duplicated_probe_count} rows with duplicated probe IDs")

    non_numeric_columns: list[str] = [
        column
        for column in group_columns
        if not atlas_df.schema[column].is_numeric()
    ]

    if non_numeric_columns:
        raise TypeError(f"Expression columns must be numeric: {non_numeric_columns}")

    null_counts: dict[str, int] = (
        atlas_df
        .select(group_columns)
        .null_count()
        .row(0, named=True)
    )

    columns_with_nulls: dict[str, int] = {
        column: count
        for column, count in null_counts.items()
        if count > 0
    }

    if columns_with_nulls:
        raise ValueError(f"Null expression values detected: {columns_with_nulls}")

    probe_ids: list[str] = (
        atlas_df
        .get_column(probe_id_column)
        .cast(pl.String)
        .to_list()
    )

    expression: np.ndarray = (
        atlas_df
        .select(group_columns)
        .to_numpy()
        .astype(np.float64)
    )

    if not np.isfinite(expression).all():
        raise ValueError("Expression matrix contains non-finite values")

    if expression_scale == "log2":
        expression = np.power(
            2.0,
            expression,
        )

    if (expression < 0).any():
        raise ValueError("Linear expression values cannot be negative")

    group_names = np.asarray(group_columns, dtype=str)

    probe_calls: list[dict[str, object]] = []

    for probe_id, probe_expression in zip(probe_ids, expression, strict=True):
        # Primary key: decreasing expression.
        # Secondary key: alphabetical group name for stable ties.
        ranking = np.lexsort(
            (
                group_names,
                -probe_expression,
            )
        )

        ranked_expression = probe_expression[ranking]
        ranked_groups = group_names[ranking]

        qualifying_boundaries: list[
            dict[str, object]
        ] = []

        for j in range(
            1,
            maximum_selected_groups + 1,
        ):
            boundary_high = float(
                ranked_expression[j - 1]
            )

            outside_expression = float(
                ranked_expression[j]
            )

            dynamic_threshold = (
                outside_expression**2 / 50.0
                + outside_expression
                + 50.0
            )

            qualifies = (
                boundary_high >= minimum_expression
                and boundary_high > dynamic_threshold
            )

            if not qualifies:
                continue

            qualifying_boundaries.append(
                {
                    "j": j,
                    "selected_groups": ranked_groups[:j].tolist(),
                    "boundary_high": boundary_high,
                    "outside_expression": outside_expression,
                    "dynamic_threshold": dynamic_threshold,
                    "dynamic_margin": boundary_high - dynamic_threshold,
                }
            )

        if not qualifying_boundaries:
            probe_calls.append(
                {
                    probe_id_column: probe_id,
                    "is_tissue_restricted": False,
                    "primary_j": 0,
                    "selected_groups": "",
                    "all_qualifying_j": "",
                    "boundary_high": None,
                    "outside_expression": None,
                    "dynamic_threshold": None,
                    "dynamic_margin": None,
                }
            )

            continue

        if selection_rule == "largest_j":
            primary_boundary = max(
                qualifying_boundaries,
                key=lambda result: int(result["j"]),
            )
        else:
            primary_boundary = min(
                qualifying_boundaries,
                key=lambda result: int(result["j"]),
            )

        selected_groups = primary_boundary[
            "selected_groups"
        ]

        probe_calls.append(
            {
                probe_id_column: probe_id,
                "is_tissue_restricted": True,
                "primary_j": primary_boundary["j"],
                "selected_groups": ";".join(
                    str(group)
                    for group in selected_groups
                ),
                "all_qualifying_j": ";".join(
                    str(result["j"])
                    for result in qualifying_boundaries
                ),
                "boundary_high": primary_boundary[
                    "boundary_high"
                ],
                "outside_expression": primary_boundary[
                    "outside_expression"
                ],
                "dynamic_threshold": primary_boundary[
                    "dynamic_threshold"
                ],
                "dynamic_margin": primary_boundary[
                    "dynamic_margin"
                ],
            }
        )

    probe_calls_df = pl.DataFrame(
        probe_calls,
        schema={
            probe_id_column: pl.String,
            "is_tissue_restricted": pl.Boolean,
            "primary_j": pl.Int64,
            "selected_groups": pl.String,
            "all_qualifying_j": pl.String,
            "boundary_high": pl.Float64,
            "outside_expression": pl.Float64,
            "dynamic_threshold": pl.Float64,
            "dynamic_margin": pl.Float64,
        },
        strict=False,
    )

    metadata_df: pl.DataFrame = atlas_df.select(probe_id_column, *metadata_columns)

    return metadata_df.join(
        probe_calls_df,
        on=probe_id_column,
        how="left",
        validate="1:1",
    )