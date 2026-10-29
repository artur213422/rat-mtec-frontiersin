from pathlib import Path

import polars as pl

from src.common.dataset_catalog.internals.types import DatasetMetadata


def prepare_isg_gene_sets(input_file: Path) -> dict[str, list[str]]:
    isgs_df: pl.DataFrame = pl.read_csv(input_file, schema_overrides={"ortholog_id": pl.String})

    required_columns: set[str] = {"target_condition", "ortholog_id"}
    missing_columns: set[str] = required_columns - set(isgs_df.columns)
    if missing_columns:
        raise ValueError(f"{input_file}: missing required columns {sorted(missing_columns)}")

    isg_members: pl.DataFrame = (
        isgs_df
        .select(
            pl.col("target_condition").cast(pl.String).str.strip_chars(),
            pl.col("ortholog_id").cast(pl.String).str.split("///"),
        )
        .explode("ortholog_id")
        .with_columns(pl.col("ortholog_id").str.strip_chars())
        .drop_nulls()
        .filter(
            (pl.col("target_condition") != "")
            & (pl.col("ortholog_id") != "")
        )
        .unique()
    )

    gene_sets: dict[str, list[str]] = {
        condition: genes
        for condition, genes in (
            isg_members
            .group_by("target_condition")
            .agg(pl.col("ortholog_id").sort())
            .sort("target_condition")
            .iter_rows()
        )
    }

    if not gene_sets:
        raise ValueError(f"No ISG gene sets found in {input_file}")

    return gene_sets


def prepare_ranked_genes(dataset: DatasetMetadata) -> tuple[pl.DataFrame, dict[str, str]]:
    tec_df: pl.DataFrame = dataset.get_standardized_orthologs_df()

    required_columns: set[str] = {
        "gene_symbol",
        "has_clean_ortholog",
        "ortholog_id",
        "ranking_statistic",
    }
    missing_columns: set[str] = required_columns - set(tec_df.columns)
    if missing_columns:
        raise ValueError(f"{dataset.dataset_id}: missing required columns {sorted(missing_columns)}")

    ranked_source: pl.DataFrame = tec_df.filter(pl.col("has_clean_ortholog"))

    rnk_df: pl.DataFrame = (
        ranked_source
        .select(
            pl.col("ortholog_id").cast(pl.String).str.strip_chars(),
            pl.col("ranking_statistic").cast(pl.Float64),
        )
        .drop_nulls()
        .filter(
            (pl.col("ortholog_id") != "")
            & pl.col("ranking_statistic").is_finite()
        )
        .sort("ranking_statistic", descending=True)
    )

    duplicate_ids: pl.DataFrame = (
        rnk_df
        .group_by("ortholog_id")
        .len()
        .filter(pl.col("len") > 1)
    )
    if duplicate_ids.height:
        raise ValueError(
            f"{dataset.dataset_id}: multiple ranked rows map to the same ortholog "
            "ID. Resolve target-side mappings or select an explicit score-collapse "
            f"rule before GSEA. Examples:\n{duplicate_ids.head(20)}"
        )

    symbol_mapping: pl.DataFrame = (
        ranked_source
        .select(
            pl.col("ortholog_id")
            .cast(pl.String)
            .str.strip_chars()
            .str.to_lowercase(),
            pl.col("gene_symbol").cast(pl.String).str.strip_chars(),
        )
        .drop_nulls()
        .filter(
            (pl.col("ortholog_id") != "")
            & (pl.col("gene_symbol") != "")
        )
        .unique()
        .group_by("ortholog_id")
        .agg(pl.col("gene_symbol").sort().str.join("///"))
    )

    return rnk_df, dict(symbol_mapping.iter_rows())


def attach_gene_symbols(output: Path, gene_symbols_by_dataset: dict[str, dict[str, str]],) -> None:
    gsea_df: pl.DataFrame = pl.read_csv(
        output,
        schema_overrides={
            "dataset_id": pl.String,
            "leading_edge_genes": pl.String,
        },
    )

    def map_symbols(row: dict[str, object]) -> str:
        dataset_id = row["dataset_id"]
        ortholog_ids = row["leading_edge_genes"]

        if not isinstance(dataset_id, str) or not isinstance(ortholog_ids, str):
            return ""

        id_to_symbol: dict[str, str] = gene_symbols_by_dataset.get(dataset_id, {})
        return ";".join(
            id_to_symbol.get(ortholog_id.lower(), ortholog_id)
            for value in ortholog_ids.split(";")
            if (ortholog_id := value.strip())
        )

    gsea_df = gsea_df.with_columns(
        pl.struct("dataset_id", "leading_edge_genes")
        .map_elements(map_symbols, return_dtype=pl.String)
        .alias("leading_edge_gene_symbols")
    )
    gsea_df.write_csv(output)


