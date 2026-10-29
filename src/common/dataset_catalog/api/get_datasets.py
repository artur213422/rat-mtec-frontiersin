from pathlib import Path

from src.common.dataset_catalog.internals.dataset_metada_parser import get_datasets
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.types.enums import CellPopulation, Species
import polars as pl

def get_mouse_tec_high() -> list[DatasetMetadata]:
    return get_datasets(Species.MOUSE, CellPopulation.MTEC_HIGH)

def get_rat_tec_high() -> DatasetMetadata:
    datasets: list[DatasetMetadata] = get_datasets(Species.RAT, CellPopulation.MTEC_HIGH)

    if len(datasets) != 1:
        raise RuntimeError("More than one dataset in get_rat_tec_high")
    return datasets[0]




def get_mouse_tec_low() -> list[DatasetMetadata]:
    return get_datasets(Species.MOUSE, CellPopulation.MTEC_LOW)


def get_rat_tec_low() -> DatasetMetadata:
    datasets: list[DatasetMetadata] = get_datasets(Species.RAT, CellPopulation.MTEC_LOW)

    if len(datasets) != 1:
        raise RuntimeError("More than one dataset in get_rat_tec_low")

    return datasets[0]

def get_mtec_low() -> list[DatasetMetadata]:
    return [
        *get_mouse_tec_low(),
        get_rat_tec_low(),
    ]


def get_mtec_high() -> list[DatasetMetadata]:
    return [
        *get_mouse_tec_high(),
        get_rat_tec_high(),
    ]

def read_file_with_ortholog_id(file: Path) -> pl.DataFrame:
    return pl.read_csv(
        file,
        schema_overrides={
            "ortholog_id": pl.String,
        },
    )