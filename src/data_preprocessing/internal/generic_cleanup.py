from polars import DataFrame


def standardize_sample_group_col_names(df: DataFrame, ctrl_cols: list[str], ko_cols: list[str]) -> DataFrame:
    missing_cols: list[str] = [
        col
        for col in [*ctrl_cols, *ko_cols]
        if col not in df.columns
    ]

    if missing_cols:
        raise ValueError(
            "Missing sample columns:\n"
            f"{missing_cols}\n\n"
            "Available columns:\n"
            f"{df.columns}"
        )

    rename_map: dict[str, str] = {}

    for index, col in enumerate(ctrl_cols, start=1):
        rename_map[col] = f"ctrl{index}"

    for index, col in enumerate(ko_cols, start=1):
        rename_map[col] = f"ko{index}"

    keep_cols: list[str] = [
        col
        for col in ["ID_REF", "gene_symbol", *ctrl_cols, *ko_cols]
        if col in df.columns
    ]

    return (
        df
        .select(keep_cols)
        .rename(rename_map)
    )