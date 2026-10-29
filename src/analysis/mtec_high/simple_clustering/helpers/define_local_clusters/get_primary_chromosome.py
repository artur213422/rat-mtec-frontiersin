import polars as pl
from src.analysis.mtec_high.simple_clustering.helpers.match_clusters.clean_orthologs import clean_orthologs

from src.common.dataset_catalog.internals.types import DatasetMetadata


_CHROMOSOME_COLUMN: str = "chromosome"
_PRIMARY_CHROMOSOME_PATTERN: str = r"^chr(?:\d+|X|Y)$"

def get_primary_chromosome_df(dataset: DatasetMetadata) -> pl.DataFrame:
    tec_df: pl.DataFrame = dataset.get_standardized_df()
    return _filter_primary_chromosomes(tec_df)


def get_primary_chromosome_orthologs_df(dataset: DatasetMetadata) -> pl.DataFrame:
    tec_df = clean_orthologs(dataset.get_standardized_df().with_columns(pl.lit(dataset.dataset_id).alias("dataset_id")))
    return _filter_primary_chromosomes(tec_df)


def _filter_primary_chromosomes(df: pl.DataFrame) -> pl.DataFrame:
    invalid_df: pl.DataFrame = df.filter(pl.col(_CHROMOSOME_COLUMN).is_not_null() & ~pl.col(_CHROMOSOME_COLUMN).str.starts_with("chr"))

    if invalid_df.height > 0:
        invalid_chromosomes: list[str] = invalid_df.get_column(_CHROMOSOME_COLUMN).unique().to_list()

        raise ValueError(f"Unexpected chromosome names: {invalid_chromosomes}")

    filtered_df: pl.DataFrame = df.filter(pl.col(_CHROMOSOME_COLUMN).str.contains(_PRIMARY_CHROMOSOME_PATTERN).fill_null(False))

    dropped_count: int = df.height - filtered_df.height

    print(f"Dropped {dropped_count}/{df.height} genes outside primary chromosomes")

    return filtered_df