import numpy as np
import numpy.typing as npt

from src.common.types.dataset import FloatMatrix
from src.exploration.general_statistics.utils.types import ReplicateStatistics



def calc_replicate_statistics(ctrl_data: FloatMatrix, ko_data: FloatMatrix,) -> ReplicateStatistics:
    mean_ctrl_correlation, minimum_ctrl_correlation = (
        _calc_pairwise_sample_correlations(
            matrix=ctrl_data,
            group_name="CTRL",
        )
    )

    mean_ko_correlation, minimum_ko_correlation = (
        _calc_pairwise_sample_correlations(
            matrix=ko_data,
            group_name="KO",
        )
    )

    return ReplicateStatistics(
        mean_ctrl_correlation=mean_ctrl_correlation,
        minimum_ctrl_correlation=minimum_ctrl_correlation,
        mean_ko_correlation=mean_ko_correlation,
        minimum_ko_correlation=minimum_ko_correlation,
    )


def _calc_pairwise_sample_correlations(matrix: FloatMatrix,group_name: str) -> tuple[float, float]:
    """
    Calculate mean and minimum Pearson correlation between samples.

    Matrix shape:
        rows = genes
        columns = replicate samples
    """
    if matrix.ndim != 2:
        raise ValueError(f"{group_name} matrix must be two-dimensional, got {matrix.ndim} dimensions.")

    gene_count: int
    sample_count: int
    gene_count, sample_count = matrix.shape

    if gene_count < 2:
        raise ValueError(f"{group_name} matrix must contain at least two genes.")

    if sample_count < 2:
        raise ValueError(f"{group_name} requires at least two replicate samples to calculate correlation.")

    if not np.isfinite(matrix).all():
        raise ValueError(f"{group_name} matrix contains NaN or infinite values.")

    sample_sds: npt.NDArray[np.float64] = np.std(
        matrix,
        axis=0,
        ddof=1,
    )

    constant_sample_indices: npt.NDArray[np.int_] = np.flatnonzero(
        sample_sds == 0.0
    )

    if constant_sample_indices.size > 0:
        raise ValueError(f"{group_name} contains constant sample columns at indices {constant_sample_indices.tolist()}. Pearson correlation is undefined for a constant sample.")

    correlation_matrix: npt.NDArray[np.float64] = np.corrcoef(
        matrix,
        rowvar=False,
    )

    upper_triangle_indices: tuple[
        npt.NDArray[np.int_],
        npt.NDArray[np.int_],
    ] = np.triu_indices(
        sample_count,
        k=1,
    )

    pairwise_correlations: npt.NDArray[np.float64] = correlation_matrix[upper_triangle_indices]


    if not np.isfinite(pairwise_correlations).all():
        raise ValueError(f"{group_name} produced undefined replicate correlations.")

    return (
        float(np.mean(pairwise_correlations)),
        float(np.min(pairwise_correlations)),
    )