import polars as pl

from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.feratures.gsea.rnk import get_rnk_from_precalulated_source
from src.common.types.enums import Species



def make_tissue_gene_sets(tsa_df: pl.DataFrame) -> dict[str, list[str]]:
    """Build strict (single-tissue) TSA gene sets for GSEA."""
    required_columns = {"gene_symbol", "primary_j", "selected_groups"}
    missing_columns = required_columns - set(tsa_df.columns)

    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    # A strict TSA is assigned to exactly one tissue group.  The input contains
    # one copy of each TSA gene for every dataset, so select only TSA-definition
    # columns and remove those repeated copies before creating the gene sets.
    tissue_members = (
        tsa_df
        .filter(pl.col("primary_j") <= 3)
        .select("gene_symbol", "selected_groups")
        .filter(
            pl.col("gene_symbol").is_not_null()
            & (pl.col("gene_symbol").str.strip_chars() != "")
            & pl.col("selected_groups").is_not_null()
            & (pl.col("selected_groups").str.strip_chars() != "")
        )
        .unique()
        .with_columns(
            pl.col("gene_symbol").str.strip_chars(),
            pl.col("selected_groups").str.split(";").alias("tissue"),
        )
        .explode("tissue")
        .with_columns(pl.col("tissue").str.strip_chars())
        .filter(pl.col("tissue") != "")
        .select("tissue", "gene_symbol")
        .unique()
    )

    return {
        tissue: genes
        for tissue, genes in (
            tissue_members
            .group_by("tissue")
            .agg(pl.col("gene_symbol").sort())
            .sort("tissue")
            .iter_rows()
        )
    }


def get_mouse_symbol_ranking(dataset: DatasetMetadata, tec_df: pl.DataFrame) -> pl.DataFrame:
    """Return a genome-wide ranking using mouse symbols for every species."""
    if dataset.species == Species.RAT:
        tec_df = (
            tec_df
            .filter(pl.col("has_clean_ortholog"))
            .select(
                pl.col("mouse_ortholog_name").alias("gene_symbol"),
                "ranking_statistic",
            )
        )

    return get_rnk_from_precalulated_source(tec_df)
