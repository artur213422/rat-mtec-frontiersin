import re

import polars as pl
import numpy as np
import numpy.typing as npt



pattern = re.compile(r"^(ctrl|ko)\d+$")

def extract_sample_names(df: pl.DataFrame) -> list[str]:
    sample_names: list[str] = []

    for name in df.columns:
        if pattern.match(name):
            sample_names.append(name)

    return sample_names


def extract_expression_data(df: pl.DataFrame, sample_ids: list[str]) -> pl.DataFrame:
    missing_columns: set[str] = set(sample_ids) - set(df.columns)

    if missing_columns:
        raise ValueError(f"Missing sample columns: {sorted(missing_columns)}")

    return  (
        df.select(sample_ids)
        .with_columns(pl.all().cast(pl.Float64))
    )

