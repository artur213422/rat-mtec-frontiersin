from dataclasses import dataclass
from pathlib import Path
import polars as pl
from src.common.paths import find_file_in_project, resolve_dir_from_root
from src.common.types.enums import Species




@dataclass(slots=True)
class DatasetMetadata:
    cel_protocol: str
    dataset_name: str
    annotation_table_file: str
    dataset_id: str
    platform_package: str
    species: Species
    ctrl_cols: list[str]
    ko_cols: list[str]

    def _get_standardized_data_file(self, root: Path = resolve_dir_from_root("data/standardized")) -> Path:
        if not root.is_dir():
            raise FileNotFoundError(f"{root}: folder not found")


        file_name: str = f"{self.dataset_id}.csv"
        path: Path = find_file_in_project(file_name, root)

        if not path.exists():
            raise FileNotFoundError(f"{file_name}: file not found")

        return path

    def get_standardized_df(self, root: Path = resolve_dir_from_root("data/standardized")) -> pl.DataFrame:
        file: Path = self._get_standardized_data_file(root)

        return pl.read_csv(
            file,
            schema_overrides={
                "ortholog_id": pl.String,
            },
            infer_schema_length=10000,
        )

    def get_standardized_orthologs_df(self, root: Path = resolve_dir_from_root("data/standardized")) -> pl.DataFrame:
        df: pl.DataFrame = self.get_standardized_df(root)

        if "has_clean_ortholog" not in df.columns:
            raise ValueError(
                f"{self.dataset_id}: missing required column 'has_clean_ortholog'"
            )

        return df.filter(pl.col("has_clean_ortholog"))


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
        return [
            dataset
            for dataset in self.get_by_species(Species.MOUSE)
            if dataset.dataset_id != "mouse_gse14365"
        ]

    def get_rat_dataset(self) -> DatasetMetadata:
        rat_dataset: list[DatasetMetadata] = self.get_by_species(Species.RAT)

        if len(rat_dataset) > 1:
            raise RuntimeError("There should be only one rat dataset")

        return rat_dataset[0]



