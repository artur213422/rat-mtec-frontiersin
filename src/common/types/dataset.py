from dataclasses import dataclass
from typing import Any

import numpy as np
import numpy.typing as npt

FloatMatrix = npt.NDArray[np.float64]

@dataclass(slots=True)
class GeneReference:
    gene_symbol: str
    ortholog_id: str | None

    @staticmethod
    def from_row(row: dict[str, Any]) -> "GeneReference":
        raw_ortholog_id: object | None = row.get("ortholog_id")

        if raw_ortholog_id is None:
            ortholog_id: str | None = None
        elif str(raw_ortholog_id).strip().lower() in ["", "none", "null", "nan"]:
            ortholog_id = None
        else:
            ortholog_id = str(raw_ortholog_id)

        return GeneReference(
            gene_symbol=str(row["gene_symbol"]),
            ortholog_id=ortholog_id,
        )