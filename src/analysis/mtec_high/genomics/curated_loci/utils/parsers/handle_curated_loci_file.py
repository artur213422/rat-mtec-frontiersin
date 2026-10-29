import polars as pl

from src.analysis.mtec_high.genomics.curated_loci.utils.configs.config import CURATED_CLUSTER_COLS
from src.analysis.mtec_high.genomics.curated_loci.utils.configs.paths import CURATED_LOCI_FILE
from src.analysis.mtec_high.genomics.curated_loci.utils.configs.types import LocusDefinitionType
from src.common.dataset_catalog.api.get_datasets import read_file_with_ortholog_id


def overwrite_to_curated_loci_file(result: pl.DataFrame, locus_definition_type: LocusDefinitionType) -> None:
    missing_columns: set[str] = set(CURATED_CLUSTER_COLS) - set(result.columns)

    if missing_columns:
        raise ValueError(f"Result is missing required columns: {sorted(missing_columns)}")

    result = result.select(CURATED_CLUSTER_COLS)

    CURATED_LOCI_FILE.parent.mkdir(parents=True, exist_ok=True)

    if not CURATED_LOCI_FILE.exists():
        result.write_csv(CURATED_LOCI_FILE)
        return

    existing_df: pl.DataFrame = read_file_with_ortholog_id(CURATED_LOCI_FILE)

    if "locus_definition_type" not in existing_df.columns:
        raise ValueError(f"Existing curated loci file does not contain 'locus_definition_type': {CURATED_LOCI_FILE}")

    existing_df = existing_df.filter(pl.col("locus_definition_type") != locus_definition_type.value)

    combined_df: pl.DataFrame = pl.concat(
        [
            existing_df.select(CURATED_CLUSTER_COLS),
            result,
        ],
        how="vertical",
    )

    combined_df.write_csv(CURATED_LOCI_FILE)

def build_curated_locus(tec_df: pl.DataFrame, member_df: pl.DataFrame, *, locus_name: str, locus_type: LocusDefinitionType, dataset_id: str) -> pl.DataFrame:
    if member_df.is_empty():
        return pl.DataFrame(
            schema={
                column: tec_df.schema.get(column, pl.String)
                for column in CURATED_CLUSTER_COLS
            }
        )

    cluster_ids: list[int] = member_df.get_column("cluster_id").unique().to_list()

    member_genes: list[str] = member_df.get_column("gene_symbol").unique().to_list()

    locus_df: pl.DataFrame = (
        tec_df
        .filter(pl.col("cluster_id").is_in(cluster_ids))
        .with_columns(
            pl.lit(locus_name).alias("locus_name"),
            pl.lit(locus_type.value).alias("locus_definition_type"),
            pl.lit(dataset_id).alias("dataset_id"),
            pl.col("gene_symbol")
            .is_in(member_genes)
            .alias("is_locus_member"),
        )
        .select(CURATED_CLUSTER_COLS)
    )

    return locus_df