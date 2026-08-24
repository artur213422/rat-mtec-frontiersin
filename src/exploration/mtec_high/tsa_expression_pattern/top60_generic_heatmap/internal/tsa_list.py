from dataclasses import dataclass
from pathlib import Path

import polars as pl

from src.common.types.dataset import GeneReference


@dataclass
class DefinedTsa:
    tissue: str
    ranked_tsa: list[GeneReference]

    @classmethod
    def from_file(cls, tissue: str, file: Path) -> "DefinedTsa":
        data: pl.DataFrame = pl.read_csv(file)

        return cls(
            tissue=tissue,
            ranked_tsa=[
                GeneReference.from_row(row)
                for row in data.iter_rows(named=True)
                if row["has_clean_ortholog"]
            ],
        )


