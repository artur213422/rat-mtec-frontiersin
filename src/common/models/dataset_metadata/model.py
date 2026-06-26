from dataclasses import dataclass
from pathlib import Path

from src.common.paths import find_file_in_project, resolve_dir_from_root
from src.common.types.enums import Species


@dataclass(slots=True)
class DatasetMetadata:
    dataset_id: str
    file: Path
    species: Species
    ctrl_cols: list[str]
    ko_cols: list[str]

    def get_standardized_data_file(self, root: Path = resolve_dir_from_root("data/standardized")) -> Path:

        if not root.is_dir():
            raise FileNotFoundError(f"{root}: Folder not found")

        file_name: str = f"{self.dataset_id}.csv"
        path: Path = find_file_in_project(file_name, root)

        if not path.exists():
            raise FileNotFoundError(f"{file_name}: File not found")

        return path



@dataclass(slots=True)
class DatasetCatalog:
    datasets: list[DatasetMetadata]


    def get_by_species(self, species: Species) -> list[DatasetMetadata]:
        return [
            dataset
            for dataset in self.datasets
            if dataset.species == species
        ]

    def get_mouse_datasets(self) -> list[DatasetMetadata]:
        return self.get_by_species(Species.MOUSE)

    def get_rat_datasets(self) -> DatasetMetadata:
        rat_dataset: list[DatasetMetadata] = self.get_by_species(Species.RAT)

        if len(rat_dataset) > 1:
            raise RuntimeError("There should be only one rat dataset")

        return rat_dataset[0]



