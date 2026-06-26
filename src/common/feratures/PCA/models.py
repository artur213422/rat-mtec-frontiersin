from dataclasses import dataclass
from enum import Enum

import polars as pl
import numpy as np
import numpy.typing as npt



@dataclass(frozen=True)
class FeatureLoading:
    genes: list[str]
    pc_names: list[str]
    gene_by_pc_loadings: npt.NDArray[np.float64]

@dataclass(frozen=True)
class SampleCoordinates:
    sample_ids: list[str]
    pc_names: list[str]
    sample_by_pc_coordinates: npt.NDArray[np.float64]


@dataclass(frozen=True)
class ExplainedVariance:
    pc_names: list[str]
    explained_variance_ratio: list[float]




@dataclass(frozen=True)
class PcaResult:
    sample_coordinates: SampleCoordinates
    explained_variance: ExplainedVariance
    feature_loadings: FeatureLoading




class Orientation(Enum):
    UNSET = 0
    GENES_BY_SAMPLES = 1
    SAMPLES_BY_GENES = 2

@dataclass(frozen=True)
class ExpressionTable:
    expression_matrix: pl.DataFrame
    genes: list[str]
    sample_ids: list[str]
    orientation: Orientation = 0

    def flip_orientation(self) -> "ExpressionTable":
        if self.orientation == Orientation.UNSET:
            raise ValueError("Cannot flip orientation because orientation is unset.")

        transposed_matrix = self.expression_matrix.transpose(
            include_header=False,
        )

        new_columns: list[str] = []
        new_orientation: Orientation = Orientation.UNSET

        if self.orientation == Orientation.GENES_BY_SAMPLES:
            new_orientation = Orientation.SAMPLES_BY_GENES
            new_columns = self.genes

        if self.orientation == Orientation.SAMPLES_BY_GENES:
            new_orientation = Orientation.GENES_BY_SAMPLES
            new_columns = self.sample_ids

        transposed_matrix.columns = new_columns

        return ExpressionTable(
            expression_matrix=transposed_matrix,
            genes=self.genes,
            sample_ids=self.sample_ids,
            orientation=new_orientation,
        )

    def to_numpy(self) -> npt.NDArray[np.float64]:
        return self.expression_matrix.to_numpy()

    def validate(self) -> None:
        if self.orientation == Orientation.UNSET:
            raise ValueError("orientation is unset")

        if self.expression_matrix.height == 0:
            raise ValueError("expression_matrix has zero rows.")

        if self.expression_matrix.width == 0:
            raise ValueError("expression_matrix has zero columns.")

        if len(self.genes) == 0:
            raise ValueError("genes list is empty.")

        if len(self.sample_ids) == 0:
            raise ValueError("sample_ids list is empty.")

        if len(set(self.genes)) != len(self.genes):
            raise ValueError("genes contains duplicate names.")

        if len(set(self.sample_ids)) != len(self.sample_ids):
            raise ValueError("sample_ids contains duplicate names.")

        if len(set(self.expression_matrix.columns)) != len(self.expression_matrix.columns):
            raise ValueError("expression_matrix contains duplicate column names.")

        allowed_dtypes = {
            pl.Float32,
            pl.Float64,
        }

        non_float_columns: list[str] = [
            column_name
            for column_name, dtype in self.expression_matrix.schema.items()
            if dtype not in allowed_dtypes
        ]

        if non_float_columns:
            raise TypeError(
                f"expression_matrix must contain only float columns. "
                f"Non-float columns found: {non_float_columns}"
            )

        null_counts = self.expression_matrix.null_count().row(0)

        columns_with_nulls: list[str] = [
            column_name
            for column_name, null_count in zip(
                self.expression_matrix.columns,
                null_counts,
                strict=True,
            )
            if null_count > 0
        ]

        if columns_with_nulls:
            raise ValueError(
                f"expression_matrix contains null values in columns: {columns_with_nulls}"
            )

        finite_checks = self.expression_matrix.select(
            [
                pl.col(column_name).is_finite().all().alias(column_name)
                for column_name in self.expression_matrix.columns
            ]
        ).row(0)

        columns_with_nan_or_inf: list[str] = [
            column_name
            for column_name, is_finite in zip(
                self.expression_matrix.columns,
                finite_checks,
                strict=True,
            )
            if not is_finite
        ]

        if columns_with_nan_or_inf:
            raise ValueError(
                f"expression_matrix contains NaN or infinite values in columns: "
                f"{columns_with_nan_or_inf}"
            )

        if self.orientation == Orientation.GENES_BY_SAMPLES:
            if len(self.genes) != self.expression_matrix.height:
                raise ValueError(
                    f"GENES_BY_SAMPLES requires one row per gene: "
                    f"{len(self.genes)} genes != {self.expression_matrix.height} rows"
                )

            if len(self.sample_ids) != self.expression_matrix.width:
                raise ValueError(
                    f"GENES_BY_SAMPLES requires one column per sample: "
                    f"{len(self.sample_ids)} samples != {self.expression_matrix.width} columns"
                )

            if self.expression_matrix.columns != self.sample_ids:
                raise ValueError(
                    "For GENES_BY_SAMPLES, matrix columns must match sample_ids.\n"
                    f"sample_ids={self.sample_ids}\n"
                    f"columns={self.expression_matrix.columns}"
                )

        elif self.orientation == Orientation.SAMPLES_BY_GENES:
            if len(self.sample_ids) != self.expression_matrix.height:
                raise ValueError(
                    f"SAMPLES_BY_GENES requires one row per sample: "
                    f"{len(self.sample_ids)} samples != {self.expression_matrix.height} rows"
                )

            if len(self.genes) != self.expression_matrix.width:
                raise ValueError(
                    f"SAMPLES_BY_GENES requires one column per gene: "
                    f"{len(self.genes)} genes != {self.expression_matrix.width} columns"
                )

            if self.expression_matrix.columns != self.genes:
                raise ValueError(
                    "For SAMPLES_BY_GENES, matrix columns must match genes.\n"
                    f"genes={self.genes[:10]}...\n"
                    f"columns={self.expression_matrix.columns[:10]}..."
                )

        else:
            raise ValueError(f"Unknown orientation: {self.orientation}")

