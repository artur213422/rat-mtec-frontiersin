from dataclasses import dataclass
from pathlib import Path
import polars as pl
from src.common.paths import resolve_dir_from_root, find_file_in_project
from src.common.types.enums import Species, CellPopulation

@dataclass(frozen=True, slots=True)
class DatasetMetadata:
    annotation_table_file: str
    dataset_name: str
    dataset_id: str
    platform_package: str
    cel_protocol: str
    species: Species
    ctrl_cols: list[str]
    ko_cols: list[str]
    cell_population: CellPopulation

    def _get_standardized_data_file(self, root: Path|None) -> Path:
        if root is None:
            if self.cell_population == CellPopulation.MTEC_HIGH:
                root = resolve_dir_from_root("data/mtec_high/standardized")
            elif self.cell_population == CellPopulation.MTEC_LOW:
                root = resolve_dir_from_root("data/mtec_low/standardized")

        file_name: str = f"{self.dataset_id}.csv"
        path: Path = find_file_in_project(file_name, root)

        if not path.exists():
            raise FileNotFoundError(f"{file_name}: file not found")

        return path

    def get_standardized_df(self, root: Path|None = None) -> pl.DataFrame:
        file: Path = self._get_standardized_data_file(root)

        return pl.read_csv(
            file,
            schema_overrides={
                "ortholog_id": pl.String,
            },
            infer_schema_length=10000,
        )

    def get_standardized_orthologs_df(self, root: Path|None = None) -> pl.DataFrame:
        df: pl.DataFrame = self.get_standardized_df(root)

        if "has_clean_ortholog" not in df.columns:
            raise ValueError(
                f"{self.dataset_id}: missing required column 'has_clean_ortholog'"
            )

        return df.filter(pl.col("has_clean_ortholog"))



