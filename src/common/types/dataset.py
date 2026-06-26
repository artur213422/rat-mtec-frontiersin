from dataclasses import dataclass
from typing import Any



@dataclass(slots=True)
class GeneReference:
    gene_symbol: str
    ortholog_id: int

    @staticmethod
    def from_row(row: dict[str, Any]) -> "GeneReference":
        return GeneReference(
            gene_symbol=str(row["gene_symbol"]),
            ortholog_id=int(row["ortholog_id"]),
        )
        