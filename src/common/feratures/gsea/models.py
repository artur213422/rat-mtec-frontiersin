from dataclasses import dataclass
from typing import Any

import polars as pl
from polars.datatypes import Float32, Float64, String


@dataclass(slots=True)
class GseaInput: # exactly two columns: gene, score
    rnk: pl.DataFrame
    gene_sets_dict: dict[str, list[str]]
    min_size: int = 7
    max_size: int = 500
    permutation_num: int = 1000

    def __post_init__(self) -> None:
        self._lowercase()
        self._validate()

    def _lowercase(self) -> None:
        gene_col: str = self.rnk.columns[0]

        self.rnk = self.rnk.with_columns(
            pl.col(gene_col).str.to_lowercase()
        )

        self.gene_sets_dict = {
            pathway: [gene.lower() for gene in genes]
            for pathway, genes in self.gene_sets_dict.items()
        }

    def _validate(self) -> None:
        if self.rnk.width != 2:
            raise ValueError(
                f"Expected rnk to have exactly 2 columns, got {self.rnk.width}."
            )

        gene_col, score_col = self.rnk.columns

        gene_dtype = self.rnk.schema[gene_col]
        score_dtype = self.rnk.schema[score_col]

        if gene_dtype != String:
            raise TypeError(
                f"First column '{gene_col}' must contain gene names (String), "
                f"got {gene_dtype}."
            )

        if score_dtype not in (Float32, Float64):
            raise TypeError(
                f"Second column '{score_col}' must contain floating-point scores, "
                f"got {score_dtype}."
            )




@dataclass(slots=True)
class GseaResult:
    term: str
    es: float
    nes: float
    nominal_p_value: float
    fdr_q_value: float
    leading_edge_genes: list[str]

@dataclass(slots=True)
class AnalysedGseaResult:
    gsea_result: list[GseaResult]
    dataset_id: str
