from statistics import mean, median

import polars as pl

from src.common.feratures.analysers.deg_extractor import is_deg
from src.common.feratures.gsea.rnk import get_genotype_expression_values
from src.common.feratures.math.p_val import calculate_log2_fold_change, calculate_raw_p_val
from src.common.types.enums import Genotypes
from src.exploration.general_statistics.utils.types import DegStatistics



def analyse_degs(data: pl.DataFrame, threshold: float = 1.0, p_val_threshold: float = 0.05,) -> DegStatistics:
    if data.is_empty():
        raise ValueError("Cannot analyse DEGs in an empty DataFrame.")

    upregulated_log2_fc: list[float] = []
    downregulated_log2_fc: list[float] = []

    for row in data.iter_rows(named=True):
        ctrl_values: list[float] = (
            get_genotype_expression_values(
                row,
                Genotypes.CTRL,
            )
        )

        ko_values: list[float] = (
            get_genotype_expression_values(
                row,
                Genotypes.KO,
            )
        )

        p_val: float = calculate_raw_p_val(ctrl_values, ko_values)

        log2_fold_change: float = calculate_log2_fold_change(ctrl_values,ko_values)

        if not is_deg(p_val=p_val, fold_change=log2_fold_change, threshold=threshold, p_val_threshold=p_val_threshold):
            continue

        if log2_fold_change > 0.0:
            upregulated_log2_fc.append(log2_fold_change)
        elif log2_fold_change < 0.0:
            downregulated_log2_fc.append(log2_fold_change)

    return _build_deg_statistics(
        upregulated_log2_fc=upregulated_log2_fc,
        downregulated_log2_fc=downregulated_log2_fc,
        total_gene_count=data.height,
    )


def _build_deg_statistics(upregulated_log2_fc: list[float], downregulated_log2_fc: list[float],total_gene_count: int) -> DegStatistics:
    if total_gene_count <= 0:
        raise ValueError("total_gene_count must be greater than zero.")

    if not upregulated_log2_fc:
        raise ValueError("No upregulated DEGs were found, so median_up_log2_fc cannot be calculated.")

    if not downregulated_log2_fc:
        raise ValueError("No downregulated DEGs were found, so median_down_log2_fc cannot be calculated.")

    upregulated_count: int = len(upregulated_log2_fc)
    downregulated_count: int = len(downregulated_log2_fc)

    total_count: int = upregulated_count + downregulated_count


    all_log2_fc: list[float] = [
        *upregulated_log2_fc,
        *downregulated_log2_fc,
    ]

    absolute_log2_fc: list[float] = [
        abs(log2_fold_change)
        for log2_fold_change in all_log2_fc
    ]

    return DegStatistics(
        upregulated_count=upregulated_count,
        downregulated_count=downregulated_count,
        total_count=total_count,

        upregulated_fraction=(upregulated_count / total_gene_count),
        downregulated_fraction=(downregulated_count / total_gene_count),
        total_fraction=(total_count / total_gene_count),

        median_abs_log2_fc=median(absolute_log2_fc),
        mean_abs_log2_fc=mean(absolute_log2_fc),
        max_abs_log2_fc=max(absolute_log2_fc),

        median_up_log2_fc=median(upregulated_log2_fc),
        median_down_log2_fc=median(downregulated_log2_fc),
    )