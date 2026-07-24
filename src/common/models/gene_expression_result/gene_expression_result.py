from dataclasses import dataclass


@dataclass(slots=True)
class GeneExpressionResult:
    dataset_id: str
    ortholog_id: str
    gene_symbol: str

    log2_ctrl_mean: float
    log2_ko_mean: float
    log2_fold_change: float

    ctrl_sd: float
    ko_sd: float

    ctrl_n: int
    ko_n: int


