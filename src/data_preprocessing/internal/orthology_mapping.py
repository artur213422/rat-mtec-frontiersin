from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import polars as pl


class Status(Enum):
    OK = "OK"
    ORTHOLOGY_ID_NOT_FOUND = "ORTHOLOGY_ID_NOT_FOUND"


@dataclass(slots=True)
class OrthologRequest:
    rat_symbols: list[str]
    mouse_symbols: list[str]
    ortholog_id: str | None

    status: Status
    message: str

class OrthologSpecies(Enum):
    MOUSE = "mouse, laboratory"
    RAT = "rat"

def split_ortholog_ids(ortholog_id: str) -> list[str]:
    return [
        part.strip()
        for part in ortholog_id.split("///")
        if part.strip() != ""
    ]

def has_orthologs(ortholog_id: str, ortholog_map: pl.DataFrame) -> bool:
    ortholog_ids: list[str] = split_ortholog_ids(ortholog_id)

    group: pl.DataFrame = ortholog_map.filter(
        pl.col("DB Class Key").cast(pl.Utf8).is_in(ortholog_ids)
    )

    has_rat: bool = group.filter(
        pl.col("Common Organism Name") == OrthologSpecies.RAT.value
    ).height > 0

    has_mouse: bool = group.filter(
        pl.col("Common Organism Name") == OrthologSpecies.MOUSE.value
    ).height > 0

    return has_rat and has_mouse

def find_ortholog_id(query_symbol: str, species: OrthologSpecies, ortholog_map: pl.DataFrame) -> str | None:
    result: pl.DataFrame = (
        ortholog_map
        .filter(
            (pl.col("Common Organism Name") == species.value)
            & (pl.col("Symbol").str.to_lowercase() == query_symbol.lower())
        )
        .select("DB Class Key")
        .unique()
    )

    if result.height == 0:
        return None

    ortholog_ids: list[str] = sorted(
        str(value)
        for value in result.to_series().to_list()
    )

    return "///".join(ortholog_ids)

def map_genes(ortholog_id: str, species: OrthologSpecies, ortholog_map: pl.DataFrame) -> list[str]:
    ortholog_ids: list[str] = split_ortholog_ids(ortholog_id)

    return (
        ortholog_map
        .filter(
            pl.col("DB Class Key").cast(pl.Utf8).is_in(ortholog_ids)
            & (pl.col("Common Organism Name") == species.value)
        )
        .select("Symbol")
        .unique()
        .to_series()
        .to_list()
    )

def find_ortholog(query_symbol: str, species: OrthologSpecies, ortholog_map: pl.DataFrame) -> OrthologRequest:
    ortholog_id: str | None  = find_ortholog_id(query_symbol, species, ortholog_map)


    if ortholog_id is None:
        return OrthologRequest(
            rat_symbols= [],
            mouse_symbols= [],
            ortholog_id=None,
            status=Status.ORTHOLOGY_ID_NOT_FOUND,
            message=f"No orthology group found for {species.value} symbol {query_symbol}",
        )

    if not has_orthologs(ortholog_id, ortholog_map):
        return OrthologRequest(
            rat_symbols= [],
            mouse_symbols= [],
            ortholog_id=None,
            status=Status.ORTHOLOGY_ID_NOT_FOUND,
            message=f"{query_symbol} orthology group does not contain both rat and mouse symbols"
        )



    rat_symbol: list[str] = map_genes(ortholog_id, OrthologSpecies.RAT, ortholog_map)
    mouse_symbol: list[str] = map_genes(ortholog_id, OrthologSpecies.MOUSE, ortholog_map)

    return OrthologRequest(
        rat_symbols=rat_symbol,
        mouse_symbols=mouse_symbol,
        ortholog_id=ortholog_id,
        status=Status.OK,
        message="",
    )


def drop_genes_with_multi_orthologs_and_save(data: pl.DataFrame, output_file: Path) -> pl.DataFrame:
    multi_ortholog_genes: pl.DataFrame = data.filter(
        pl.col("ortholog_id")
        .fill_null("")
        .str.contains("///", literal=True)
    )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    multi_ortholog_genes.write_csv(output_file)

    cleaned_data: pl.DataFrame = data.filter(
        ~pl.col("ortholog_id")
        .fill_null("")
        .str.contains("///", literal=True)
    )

    return cleaned_data