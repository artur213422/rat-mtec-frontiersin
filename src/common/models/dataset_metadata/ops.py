import tomllib
from pathlib import Path
from typing import Any

from src.common.models.dataset_metadata.model import DatasetMetadata, DatasetCatalog
from src.common.paths import find_file_in_project
from src.common.types.enums import Species


DATASET_CONFIG_FILE_NAME: str = "datasets.toml"


def build_dataset_catalog() -> DatasetCatalog:
    dataset_config_path: Path = find_file_in_project(DATASET_CONFIG_FILE_NAME)

    with dataset_config_path.open("rb") as file:
        raw_config: dict[str, Any] = tomllib.load(file)

    raw_datasets: list[dict[str, Any]] = raw_config["dataset"]

    datasets: list[DatasetMetadata] = [
        _parse_dataset_metadata(raw_dataset)
        for raw_dataset in raw_datasets
    ]

    return DatasetCatalog(datasets=datasets)


def _parse_dataset_metadata(raw_dataset: dict[str, Any]) -> DatasetMetadata:
    return DatasetMetadata(
        dataset_id =raw_dataset["dataset_id"],
        file=Path(str(raw_dataset["file"])),
        species=Species(str(raw_dataset["species"])),
        ctrl_cols=[str(col) for col in raw_dataset["control_samples"]],
        ko_cols=[str(col) for col in raw_dataset["ko_samples"]],
    )