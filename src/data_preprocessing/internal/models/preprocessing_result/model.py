from math import inf
from dataclasses import dataclass


@dataclass(slots=True)
class PreprocessingResult:
    dataset_name: str
    total_gene_count: int = inf
    dropped_invalid_count: int = inf
    dropped_duplicates_count: int = inf
    unresolved_identifiers: int = inf
    cleaned_genes_count: int = inf

