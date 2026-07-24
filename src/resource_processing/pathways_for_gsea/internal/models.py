
from dataclasses import dataclass, field

from src.reference_processing.internal.enums import Namespace


@dataclass
class GoPathway:
    namespace: Namespace  # biological_process / molecular_function / cellular_component
    go_term: str        # GO:0019882
    name: str
    genes: list[str] = field(default_factory=list)


