from typing import Any

import polars as pl

from src.common.dataframe_ops import get_genotype_expression_values
from src.common.feratures.math.p_val import calculate_log2_fold_change, calculate_raw_p_val

from src.common.types.enums import Genotypes



def genotype_comparison(data: pl.DataFrame) -> pl.DataFrame:
    rows: list[dict[str, Any]] = []

    for row in data.iter_rows(named=True):
        ctrl_values: list[float] = get_genotype_expression_values(row, Genotypes.CTRL)

        ko_values: list[float] = get_genotype_expression_values(row, Genotypes.KO)

        fc: float = calculate_log2_fold_change(
            ctrl_values,
            ko_values,
        )

        p_value: float = calculate_raw_p_val(ctrl_values,ko_values)

        row["log2_fold_change"] = fc
        row["p_value"] = p_value

        rows.append(row)

    return pl.DataFrame(rows)



