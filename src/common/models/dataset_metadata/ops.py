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
        if raw_dataset["dataset_id"] != "mouse_gse14365"
    ]

    return DatasetCatalog(datasets=datasets)


def _parse_dataset_metadata(raw_dataset: dict[str, Any]) -> DatasetMetadata:
    return DatasetMetadata(
        cel_protocol=raw_dataset["cel_protocol"],
        dataset_name=raw_dataset["dataset_name"],
        annotation_table_file=raw_dataset["annotation_table"],
        dataset_id =raw_dataset["dataset_id"],
        platform_package=raw_dataset["platform_package"],
        species=Species(str(raw_dataset["species"])),
        ctrl_cols=[str(col) for col in raw_dataset["control_samples"]],
        ko_cols=[str(col) for col in raw_dataset["ko_samples"]],
    )