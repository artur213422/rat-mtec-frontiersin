import numpy as np
import numpy.typing as npt

from src.common.types.dataset import FloatMatrix
from src.exploration.general_statistics.utils.types import ExpressionStatistics



def calc_expression_statistics(matrix: FloatMatrix) -> ExpressionStatistics:
    if matrix.ndim != 2:
        raise ValueError(f"Expression matrix must be two-dimensional, got {matrix.ndim} dimensions.")

    gene_count: int
    sample_count: int
    gene_count, sample_count = matrix.shape

    if gene_count == 0:
        raise ValueError("Expression matrix contains no genes.")

    if sample_count < 2:
        raise ValueError("At least two samples are required to calculate gene-level SD.")

    if not np.isfinite(matrix).all():
        raise ValueError("Expression matrix contains NaN or infinite values.")

    gene_means: FloatMatrix = np.mean(
        matrix,
        axis=1,
        dtype=np.float64,
    )

    gene_sds: FloatMatrix = np.std(
        matrix,
        axis=1,
        ddof=1,
        dtype=np.float64,
    )

    # CV is only meaningful for genes with positive mean expression.
    valid_cv_mask: npt.NDArray[np.bool_] = gene_means > 0.0

    if not np.any(valid_cv_mask):
        raise ValueError("Cannot calculate gene CV because no gene has positive mean expression.")

    gene_cvs: FloatMatrix = gene_sds[valid_cv_mask] / gene_means[valid_cv_mask]


    return ExpressionStatistics(
        mean_expression=float(np.mean(matrix)),
        median_expression=float(np.median(matrix)),

        mean_gene_sd=float(np.mean(gene_sds)),
        median_gene_sd=float(np.median(gene_sds)),

        mean_gene_cv=float(np.mean(gene_cvs)),
        median_gene_cv=float(np.median(gene_cvs)),

        expression_q25=float(np.quantile(matrix, 0.25)),
        expression_q75=float(np.quantile(matrix, 0.75)),
    )