from dataclasses import dataclass


@dataclass(slots=True)
class GeneExpressionResult:
    dataset_id: str
    ortholog_id: int
    gene_symbol: str

    ctrl_mean: float
    ko_mean: float
    fold_change: float

    ctrl_sd: float
    ko_sd: float

    ctrl_n: int
    ko_n: int


