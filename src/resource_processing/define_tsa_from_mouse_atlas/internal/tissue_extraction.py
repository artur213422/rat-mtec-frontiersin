from __future__ import annotations

import json
import re
from pathlib import Path

import polars as pl


def extract_list_from_json(json_file: Path) -> list[str]:
    with json_file.open("r", encoding="utf-8") as file:
        data = json.load(file)

    if not isinstance(data, list):
        raise TypeError(f"Expected JSON list in {json_file}, got {type(data).__name__}")

    if not all(isinstance(item, str) for item in data):
        raise TypeError(f"Expected list[str] in {json_file}")

    return data




_REPLICATE_SUFFIX = re.compile(r"_\d+$")


def remove_replicate_suffix(col_name: str) -> str:
    return _REPLICATE_SUFFIX.sub("", col_name)


def find_potential_match(tissue: str, col_names: list[str]) -> list[str]:
    """
    Find matrix columns matching a tissue name.

    Example:
        tissue='liver'
        col_names=['liver_1', 'liver_2', 'kidney_1']
        -> ['liver_1', 'liver_2']
    """
    return [
        col
        for col in col_names
        if remove_replicate_suffix(col) == tissue
    ]


def get_selected_tissue_columns(
    picked_tissues: list[str],
    all_tissue_cols: list[str],
) -> list[str]:
    selected_cols: list[str] = []

    for tissue in picked_tissues:
        # Allows picked_tissues.json to contain either "liver" or exact "liver_1"
        if tissue in all_tissue_cols:
            selected_cols.append(tissue)
            continue

        matches = find_potential_match(tissue, all_tissue_cols)

        if not matches:
            raise ValueError(f"No matching tissue columns found for: {tissue}")

        selected_cols.extend(matches)

    return selected_cols


def extract_tissue_matrix(
    expression_matrix_file: Path,
    picked_tissues: list[str],
    all_tissue_cols: list[str],
) -> pl.DataFrame:
    matrix: pl.DataFrame = pl.read_csv(expression_matrix_file)

    selected_tissue_cols: list[str] = get_selected_tissue_columns(
        picked_tissues=picked_tissues,
        all_tissue_cols=all_tissue_cols,
    )

    cols_to_keep: list[str] = ["gene_symbol", *selected_tissue_cols]

    missing_cols: list[str] = [
        col for col in cols_to_keep if col not in matrix.columns
    ]

    if missing_cols:
        raise ValueError(f"Columns not found in expression matrix: {missing_cols}")

    return matrix.select(cols_to_keep)