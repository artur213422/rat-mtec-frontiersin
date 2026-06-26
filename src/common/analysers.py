from statistics import stdev, mean

import polars as pl
import re

from src.common.dataframe_ops import get_genotype_expression_values
from src.common.models.gene_expression_result.gene_expression_result import GeneExpressionResult
from src.common.feratures.gene_pattern_extraction.models import GeneReference
from src.common.types.enums import Genotypes



def genotype_comparison(data: pl.DataFrame, ortholog_id: int, dataset_id: str) -> GeneExpressionResult:
    result: pl.DataFrame = data.filter(
        pl.col("ortholog_id") == ortholog_id
    )

    if result.height != 1:
        raise ValueError(
            f"Expected exactly 1 row for ortholog_id={ortholog_id}, "
            f"found {result.height}"
        )

    row: dict[str, object] = result.row(0, named=True)

    ctrl_values: list[float] = get_genotype_expression_values(
        row,
        Genotypes.CTRL,
    )

    ko_values: list[float] = get_genotype_expression_values(
        row,
        Genotypes.KO,
    )

    ctrl_mean: float = mean(ctrl_values)
    ko_mean: float = mean(ko_values)

    return GeneExpressionResult(
        dataset_id=dataset_id,
        ortholog_id=ortholog_id,
        gene_symbol=str(row["gene_symbol"]),
        ctrl_mean=ctrl_mean,
        ko_mean=ko_mean,
        fold_change=ko_mean - ctrl_mean,
        ctrl_sd=stdev(ctrl_values) if len(ctrl_values) > 1 else 0.0,
        ko_sd=stdev(ko_values) if len(ko_values) > 1 else 0.0,
        ctrl_n=len(ctrl_values),
        ko_n=len(ko_values),
    )




