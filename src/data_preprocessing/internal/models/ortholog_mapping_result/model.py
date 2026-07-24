from math import inf
from dataclasses import dataclass


@dataclass(slots=True)
class OrthologMappingResult:
    dataset_name: str
    total_genes: int = 0
    genes_with_no_orthology: int = 0
    rat_polymorphism_count: int = 0
    mouse_polymorphism_count: int = 0
    genes_with_multiple_orthologs: int = 0
    genes_with_clean_ortholog_count: int = 0
