from statistics import fmean

import polars as pl

from src.common.dataframe_ops import get_genotype_expression_values
from src.common.types.enums import Genotypes


def drop_low_expressing_genes(tec_df: pl.DataFrame, threshold: float) -> pl.DataFrame:
    keep_rows: list[bool] = []
    dropped: int = 0

    for row in tec_df.iter_rows(named=True):
        ctrl_values: list[float] = get_genotype_expression_values(
            row=row,
            gt=Genotypes.CTRL,
        )
        ko_values: list[float] = get_genotype_expression_values(
            row=row,
            gt=Genotypes.KO,
        )

        if not ctrl_values:
            raise ValueError("No CTRL expression columns were found")

        if not ko_values:
            raise ValueError("No KO expression columns were found")

        ctrl_mean: float = fmean(ctrl_values)
        ko_mean: float = fmean(ko_values)

        is_low_expression: bool = (
                ctrl_mean <= threshold
                and ko_mean <= threshold
        )

        keep_rows.append(not is_low_expression)

        if is_low_expression:
            dropped += 1

    filtered_df: pl.DataFrame = tec_df.filter(
        pl.Series(
            name="_keep",
            values=keep_rows,
            dtype=pl.Boolean,
        )
    )

    print(
        f"Dropped {dropped:,} genes with CTRL and KO means "
        f"below {threshold}"
    )

    return filtered_df