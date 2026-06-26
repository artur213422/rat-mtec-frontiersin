from dataclasses import dataclass

from src.common.types.dataset import GeneReference

@dataclass(slots=True)
class RankedGene:
    gene_reference: GeneReference
    log2_fold_change: float
    p_val: float


@dataclass
class TopGenesResult:
    top_upregulated: list[RankedGene]
    top_downregulated: list[RankedGene]

