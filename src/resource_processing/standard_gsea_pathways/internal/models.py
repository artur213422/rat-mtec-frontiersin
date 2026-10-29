
from dataclasses import dataclass, field

from src.resource_processing.standard_gsea_pathways.internal.enums import Namespace


@dataclass
class GoPathway:
    namespace: Namespace  # biological_process / molecular_function / cellular_component
    go_term: str        # GO:0019882
    name: str
    parent_go_terms: list[str] = field(default_factory=list)
    genes: list[str] = field(default_factory=list)


