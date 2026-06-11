from dataclasses import dataclass, field
from enum import Enum

import polars as pl

from data.shared_sdk.gene_symbol_cleaning_tools import find_match

class Status(Enum):
    OK = "OK"
    ORTHOLOGY_ID_NOT_FOUND = "ORTHOLOGY_ID_NOT_FOUND"
    NO_ORTHOLOGY_FOUND = "NO_ORTHOLOGY_FOUND"


@dataclass
class OrthologResult:
    rat_symbols: list[str]
    mouse_symbols: list[str]
    ortholog_id: int | None

    status: Status
    message: str

class Species(Enum):
    MOUSE = "mouse, laboratory"
    RAT = "rat"


def has_orthologs(ortholog_id: int, ortholog_map: pl.DataFrame) -> bool:
    group = ortholog_map.filter(
        pl.col("DB Class Key") == ortholog_id
    )

    has_rat = group.filter(
        pl.col("Common Organism Name") == Species.RAT.value
    ).height > 0

    has_mouse = group.filter(
        pl.col("Common Organism Name") == Species.MOUSE.value
    ).height > 0

    return has_rat and has_mouse

def find_ortholog_id(query_symbol: str, species: Species, ortholog_map: pl.DataFrame) -> int | None:
    result = (
        ortholog_map
        .filter(
            (pl.col("Common Organism Name") == species.value) &
            (pl.col("Symbol").str.to_lowercase() == query_symbol.lower())
        )
        .select("DB Class Key")
        .limit(1)
    )

    if result.height == 0:
        return None

    return int(result.item())

def map_genes(ortholog_id: int, species: Species,  ortholog_map: pl.DataFrame) -> list[str]:
    return (
        ortholog_map
        .filter(
            (pl.col("DB Class Key") == ortholog_id) &
            (pl.col("Common Organism Name") == species.value)
        )
        .select(pl.col("Symbol").str.to_lowercase())
        .unique()
        .to_series()
        .to_list()
    )

def find_ortholog(query_symbol: str, species: Species, ortholog_map: pl.DataFrame) -> OrthologResult:
    ortholog_id: int | None  = find_ortholog_id(query_symbol, species, ortholog_map)


    if ortholog_id is None:
        return OrthologResult(
            rat_symbols= [],
            mouse_symbols= [],
            ortholog_id=None,
            status=Status.ORTHOLOGY_ID_NOT_FOUND,
            message=f"No orthology group found for {species.value} symbol {query_symbol}",
        )

    if not has_orthologs(ortholog_id, ortholog_map):
        return OrthologResult(
            rat_symbols= [],
            mouse_symbols= [],
            ortholog_id=None,
            status=Status.ORTHOLOGY_ID_NOT_FOUND,
            message=f"{query_symbol} has no rat or mouse orthology",
        )



    rat_symbol: list[str] = map_genes(ortholog_id, Species.RAT, ortholog_map)
    mouse_symbol: list[str] = map_genes(ortholog_id, Species.MOUSE, ortholog_map)

    return OrthologResult(
        rat_symbols=rat_symbol,
        mouse_symbols=mouse_symbol,
        ortholog_id=ortholog_id,
        status=Status.OK,
        message="",
    )