import polars as pl
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import numpy.typing as npt

from src.common.feratures.PCA.models import ExpressionTable, Orientation, PcaResult, SampleCoordinates, \
    ExplainedVariance, FeatureLoading
from src.common.feratures.PCA.utils import extract_expression_data, extract_sample_names


def prepare_pca_input(data: pl.DataFrame) -> ExpressionTable:
    sample_names: list[str] = extract_sample_names(data)

    if len(sample_names) < 5:
        raise ValueError(f"Got less than 5 samples. Is that valid? Sample names: {sample_names}")

    if "gene_symbol" not in data.columns:
        raise ValueError("Missing required column: gene_symbol")

    genes: list[str] = data["gene_symbol"].to_list()

    if len(genes) == 0:
        raise ValueError("Gene list is empty")

    matrix: pl.DataFrame = extract_expression_data(data, sample_names)

    return ExpressionTable(
        expression_matrix=matrix,
        genes=genes,
        sample_ids=sample_names,
        orientation=Orientation.GENES_BY_SAMPLES
    )


def run_pca_analysis(expression_table: ExpressionTable) -> PcaResult:
    ########  Validation  ########
    expression_table.validate()

    if  Orientation.SAMPLES_BY_GENES != expression_table.orientation:
        expression_table = expression_table.flip_orientation()

    expression_table.validate()

    n_samples: int = len(expression_table.sample_ids)
    n_genes: int = len(expression_table.genes)

    n_components: int = min(n_samples - 1, n_genes)

    ########  Analyse  ########
    x: npt.NDArray[np.float64] = expression_table.expression_matrix.to_numpy()

    x_scaled: npt.NDArray[np.float64] = StandardScaler().fit_transform(x)

    pca = PCA(n_components=n_components)
    coordinates: npt.NDArray[np.float64] = pca.fit_transform(x_scaled)

    pc_names: list[str] = [f"PC{i + 1}" for i in range(n_components)]

    sample_coordinates = SampleCoordinates(
        sample_ids=expression_table.sample_ids,
        pc_names=pc_names,
        sample_by_pc_coordinates=coordinates,
    )

    explained_variance = ExplainedVariance(
        pc_names=pc_names,
        explained_variance_ratio=pca.explained_variance_ratio_.tolist(),
    )

    feature_loadings = FeatureLoading(
        genes=expression_table.genes,
        pc_names=pc_names,
        gene_by_pc_loadings=pca.components_.T,
    )

    return PcaResult(
        sample_coordinates=sample_coordinates,
        explained_variance=explained_variance,
        feature_loadings=feature_loadings,
    )


def