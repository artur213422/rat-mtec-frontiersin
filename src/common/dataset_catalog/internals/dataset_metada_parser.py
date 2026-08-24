from pathlib import Path
from typing import Any, Final
import tomllib

from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.paths import find_file_in_project
from src.common.types.enums import CellPopulation, Species


MTEC_HIGH_DATASET_CONFIG_FILE_NAME: Final[str] = (
    "mtec_high_datasets.toml"
)
MTEC_LOW_DATASET_CONFIG_FILE_NAME: Final[str] = (
    "mtec_low_datasets.toml"
)

MTEC_HIGH_DATASET_CONFIG_PATH: Final[Path] = find_file_in_project(
    MTEC_HIGH_DATASET_CONFIG_FILE_NAME
)
MTEC_LOW_DATASET_CONFIG_PATH: Final[Path] = find_file_in_project(
    MTEC_LOW_DATASET_CONFIG_FILE_NAME
)

_CONFIG_FILE_BY_POPULATION: Final[dict[CellPopulation, Path]] = {
    CellPopulation.MTEC_HIGH: MTEC_HIGH_DATASET_CONFIG_PATH,
    CellPopulation.MTEC_LOW: MTEC_LOW_DATASET_CONFIG_PATH,
}


def _get_config_path(
    population: CellPopulation,
) -> Path:
    try:
        return _CONFIG_FILE_BY_POPULATION[population]
    except KeyError as error:
        raise ValueError(
            f"Unsupported cell population: {population!r}"
        ) from error


def _read_toml_datasets(
    config_path: Path,
) -> list[dict[str, Any]]:
    """Read and validate [[dataset]] entries from a TOML file."""
    with config_path.open(mode="rb") as file:
        config: dict[str, Any] = tomllib.load(file)

    raw_datasets: object = config.get("dataset")

    if not isinstance(raw_datasets, list):
        raise ValueError(
            f"{config_path} must contain one or more "
            "[[dataset]] sections."
        )

    datasets: list[dict[str, Any]] = []

    for index, raw_dataset in enumerate(raw_datasets):
        if not isinstance(raw_dataset, dict):
            raise TypeError(
                f"Dataset entry {index} in {config_path} "
                "must be a TOML table."
            )

        datasets.append(raw_dataset)

    return datasets


def _require_string(
    dataset: dict[str, Any],
    key: str,
    *,
    config_path: Path,
) -> str:
    """Read a required string field from a dataset table."""
    value: object = dataset.get(key)

    if not isinstance(value, str) or not value.strip():
        dataset_id: object = dataset.get(
            "dataset_id",
            "<unknown>",
        )

        raise ValueError(
            f"Dataset {dataset_id!r} in {config_path} "
            f"must define a non-empty string field {key!r}."
        )

    return value


def _require_string_list(
    dataset: dict[str, Any],
    key: str,
    *,
    config_path: Path,
) -> list[str]:
    """Read a required list of strings from a dataset table."""
    value: object = dataset.get(key)

    if not isinstance(value, list):
        dataset_id: object = dataset.get(
            "dataset_id",
            "<unknown>",
        )

        raise ValueError(
            f"Dataset {dataset_id!r} in {config_path} "
            f"must define {key!r} as a list."
        )

    if not all(isinstance(item, str) for item in value):
        dataset_id = dataset.get(
            "dataset_id",
            "<unknown>",
        )

        raise TypeError(
            f"All entries in {key!r} for dataset "
            f"{dataset_id!r} must be strings."
        )

    return list(value)


def _parse_species(
    value: str,
    *,
    dataset_id: str,
) -> Species:
    """Convert a TOML species string into a Species enum."""
    try:
        return Species(value.lower())
    except ValueError as error:
        supported: str = ", ".join(
            member.value for member in Species
        )

        raise ValueError(
            f"Dataset {dataset_id!r} has unsupported species "
            f"{value!r}. Supported values: {supported}."
        ) from error


def get_datasets(
    species: Species,
    population: CellPopulation,
) -> list[DatasetMetadata]:
    """
    Load datasets for one species and cell population.

    The cell population determines which TOML file is loaded.
    Entries belonging to other species are ignored.
    """
    config_path: Path = _get_config_path(population)
    raw_datasets: list[dict[str, Any]] = (
        _read_toml_datasets(config_path)
    )

    datasets: list[DatasetMetadata] = []

    for raw_dataset in raw_datasets:
        dataset_id: str = _require_string(
            raw_dataset,
            "dataset_id",
            config_path=config_path,
        )

        dataset_species: Species = _parse_species(
            _require_string(
                raw_dataset,
                "species",
                config_path=config_path,
            ),
            dataset_id=dataset_id,
        )

        if dataset_species is not species:
            continue

        dataset: DatasetMetadata = DatasetMetadata(
            annotation_table_file=_require_string(
                raw_dataset,
                "annotation_table",
                config_path=config_path,
            ),
            dataset_name=_require_string(
                raw_dataset,
                "dataset_name",
                config_path=config_path,
            ),
            dataset_id=dataset_id,
            platform_package=_require_string(
                raw_dataset,
                "platform_package",
                config_path=config_path,
            ),
            cel_protocol=_require_string(
                raw_dataset,
                "cel_protocol",
                config_path=config_path,
            ),
            species=dataset_species,
            ctrl_cols=_require_string_list(
                raw_dataset,
                "control_samples",
                config_path=config_path,
            ),
            ko_cols=_require_string_list(
                raw_dataset,
                "ko_samples",
                config_path=config_path,
            ),
            cell_population=population,
        )

        datasets.append(dataset)

    return datasets