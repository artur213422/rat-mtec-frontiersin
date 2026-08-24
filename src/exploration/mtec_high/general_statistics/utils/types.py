import math
from dataclasses import dataclass

from src.common.feratures.analysers.deg_models import TopGenesResult


def _validate_finite_float(value: float, field_name: str,) -> None:
    if not math.isfinite(value):
        raise ValueError(f"{field_name} must be finite, got {value}.")


@dataclass(frozen=True, slots=True)
class ExpressionStatistics:
    mean_expression: float
    median_expression: float

    mean_gene_sd: float
    median_gene_sd: float

    mean_gene_cv: float
    median_gene_cv: float

    expression_q25: float
    expression_q75: float

    def __post_init__(self) -> None:
        for field_name, value in (
            ("mean_expression", self.mean_expression),
            ("median_expression", self.median_expression),
            ("mean_gene_sd", self.mean_gene_sd),
            ("median_gene_sd", self.median_gene_sd),
            ("mean_gene_cv", self.mean_gene_cv),
            ("median_gene_cv", self.median_gene_cv),
            ("expression_q25", self.expression_q25),
            ("expression_q75", self.expression_q75),
        ):
            _validate_finite_float(value=value, field_name=field_name)

        if self.mean_gene_sd < 0:
            raise ValueError("mean_gene_sd cannot be negative.")

        if self.median_gene_sd < 0:
            raise ValueError("median_gene_sd cannot be negative.")

        if self.expression_q25 > self.expression_q75:
            raise ValueError("expression_q25 cannot exceed expression_q75.")


@dataclass(frozen=True, slots=True)
class DegStatistics:
    upregulated_count: int
    downregulated_count: int
    total_count: int

    upregulated_fraction: float
    downregulated_fraction: float
    total_fraction: float

    median_abs_log2_fc: float
    mean_abs_log2_fc: float
    max_abs_log2_fc: float

    median_up_log2_fc: float
    median_down_log2_fc: float

    def __post_init__(self) -> None:
        expected_total: int = (self.upregulated_count + self.downregulated_count)

        if self.total_count != expected_total:
            raise ValueError("total_count must equal upregulated_count + downregulated_count.")

        for field_name, count in (
            ("upregulated_count", self.upregulated_count),
            ("downregulated_count", self.downregulated_count),
            ("total_count", self.total_count),
        ):
            if count < 0:
                raise ValueError(f"{field_name} cannot be negative.")

        for field_name, fraction in (
            ("upregulated_fraction", self.upregulated_fraction),
            ("downregulated_fraction", self.downregulated_fraction),
            ("total_fraction", self.total_fraction),
        ):
            _validate_finite_float(value=fraction, field_name=field_name)

            if not 0.0 <= fraction <= 1.0:
                raise ValueError(f"{field_name} must be between 0 and 1.")

        for field_name, value in (
            ("median_abs_log2_fc", self.median_abs_log2_fc),
            ("mean_abs_log2_fc", self.mean_abs_log2_fc),
            ("max_abs_log2_fc", self.max_abs_log2_fc),
            ("median_up_log2_fc", self.median_up_log2_fc),
            ("median_down_log2_fc", self.median_down_log2_fc),
        ):
            _validate_finite_float(
                value=value,
                field_name=field_name,
            )

        if self.median_abs_log2_fc < 0:
            raise ValueError("median_abs_log2_fc cannot be negative." )

        if self.mean_abs_log2_fc < 0:
            raise ValueError("mean_abs_log2_fc cannot be negative.")

        if self.max_abs_log2_fc < 0:
            raise ValueError("max_abs_log2_fc cannot be negative.")

        if self.upregulated_count == 0:
            raise ValueError("median_up_log2_fc is undefined because there are no upregulated DEGs.")

        if self.downregulated_count == 0:
            raise ValueError("median_down_log2_fc is undefined because there are no downregulated DEGs.")

        if self.median_up_log2_fc <= 0:
            raise ValueError("median_up_log2_fc must be positive.")

        if self.median_down_log2_fc >= 0:
            raise ValueError("median_down_log2_fc must be negative.")


@dataclass(frozen=True, slots=True)
class ReplicateStatistics:
    mean_ctrl_correlation: float
    minimum_ctrl_correlation: float

    mean_ko_correlation: float
    minimum_ko_correlation: float

    def __post_init__(self) -> None:
        for field_name, value in (
                ("mean_ctrl_correlation", self.mean_ctrl_correlation),
                ("minimum_ctrl_correlation", self.minimum_ctrl_correlation,
            ),
                ("mean_ko_correlation", self.mean_ko_correlation),
                ("minimum_ko_correlation",self.minimum_ko_correlation,),
        ):
            _validate_finite_float(value=value, field_name=field_name)

            if not -1.0 <= value <= 1.0:
                raise ValueError(f"{field_name} must be between -1 and 1.")

        if self.minimum_ctrl_correlation > self.mean_ctrl_correlation:
            raise ValueError("minimum_ctrl_correlation cannot exceed mean_ctrl_correlation.")

        if self.minimum_ko_correlation> self.mean_ko_correlation:
            raise ValueError("minimum_ko_correlation cannot exceed mean_ko_correlation.")


@dataclass(frozen=True, slots=True)
class DatasetStatistics:
    dataset_id: str
    degs: DegStatistics
    total_unique_genes: int

    ctrl: ExpressionStatistics
    ko: ExpressionStatistics
    combined: ExpressionStatistics

    replicate_statistics: ReplicateStatistics
    top_15: TopGenesResult

    def __post_init__(self) -> None:
        if not self.dataset_id.strip():
            raise ValueError("dataset_id cannot be empty.")

        if self.total_unique_genes <= 0:
            raise ValueError("total_unique_genes must be greater than zero.")

        if self.degs.total_count > self.total_unique_genes:
            raise ValueError("DEG count cannot exceed total_unique_genes.")