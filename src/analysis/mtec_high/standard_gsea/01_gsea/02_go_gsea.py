from __future__ import annotations

from pathlib import Path

import polars as pl

from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.excel import create_excel_table
from src.common.feratures.gsea.checkpoint import run_gsea_checkpointed
from src.common.feratures.gsea.models import AnalysedGseaResult, GseaInput, GseaResult
from src.common.feratures.gsea.rnk import get_rnk_from_precalulated_source
from src.common.feratures.gsea.save_result import save_gsea_result
from src.common.file_io import save_csv_to_excel
from src.common.paths import PROJECT_ROOT

##########################################################################################
#                                      SETTINGS
##########################################################################################

GO_NAMESPACES: tuple[str, ...] = (
    "biological_process",
    "molecular_function",
    "cellular_component",
)

go_pathways_file: Path = PROJECT_ROOT / "resources/go/go_pathway_genes.csv"

output_folder: Path = PROJECT_ROOT / "data/mtec_high/derived_data/standard_gsea/both_sided_gsea"

##########################################################################################
#                                  GENE-SET HELPERS
##########################################################################################

def create_go_gene_sets_dict(go_df: pl.DataFrame, *, namespace: str, species: str) -> dict[str, list[str]]:
    required_columns: set[str] = {
        "species",
        "namespace",
        "go_term",
        "name",
        "gene_symbol",
    }

    missing_columns: set[str] = required_columns - set(go_df.columns)

    if missing_columns:
        raise ValueError(f"GO pathway table is missing columns: {sorted(missing_columns)}")

    namespace_df: pl.DataFrame = (
        go_df
        .select(
            pl.col("species")
            .cast(pl.String)
            .str.to_lowercase()
            .str.strip_chars(),

            pl.col("namespace")
            .cast(pl.String)
            .str.to_lowercase()
            .str.strip_chars(),

            pl.col("go_term")
            .cast(pl.String)
            .str.strip_chars(),

            pl.col("name")
            .cast(pl.String)
            .str.strip_chars(),

            pl.col("gene_symbol")
            .cast(pl.String)
            .str.strip_chars(),
        )
        .filter(
            (pl.col("species") == species.lower())
            & (pl.col("namespace") == namespace)
            & pl.col("go_term").is_not_null()
            & (pl.col("go_term") != "")
            & pl.col("name").is_not_null()
            & (pl.col("name") != "")
            & pl.col("gene_symbol").is_not_null()
            & (pl.col("gene_symbol") != "")
        )
        .unique(
            subset=[
                "go_term",
                "gene_symbol",
            ]
        )
    )

    if namespace_df.is_empty():
        raise ValueError(f"No GO pathways found for namespace={namespace!r}, species={species!r}.")

    gene_sets: dict[str, list[str]] = {}

    for row in (
        namespace_df
        .group_by(
            "go_term",
            "name",
        )
        .agg(
            pl.col("gene_symbol")
            .unique()
            .sort()
            .alias("genes")
        )
        .sort("go_term")
        .iter_rows(named=True)
    ):
        go_term: str = str(row["go_term"])
        pathway_name: str = str(row["name"])

        # GO ID guarantees uniqueness even when pathway names repeat.
        pathway_key: str = f"{go_term} | {pathway_name}"

        genes: list[str] = [
            str(gene)
            for gene in row["genes"]
        ]

        gene_sets[pathway_key] = genes

    return gene_sets

##########################################################################################
#                                       READ DATA
##########################################################################################

standardized_folder: Path = PROJECT_ROOT / "data/mtec_high/standardized"

datasets: list[DatasetMetadata] = get_mtec_high()

go_pathways_df: pl.DataFrame = pl.read_csv(go_pathways_file, infer_schema_length=10_000)

gene_sets_by_namespace: dict[
    tuple[str, str],
    dict[str, list[str]],
] = {
    (species, namespace): create_go_gene_sets_dict(
        go_pathways_df,
        namespace=namespace,
        species=species,
    )
    for species in {dataset.species.value for dataset in datasets}
    for namespace in GO_NAMESPACES
}

for (species, namespace), gene_sets in gene_sets_by_namespace.items():
    print(f"{species} / {namespace}: {len(gene_sets):,} gene sets")

##########################################################################################
#                                       RUN GSEA
##########################################################################################

for namespace in GO_NAMESPACES:
    analysed_results: list[AnalysedGseaResult] = []

    print()
    print(f"Running namespace: {namespace}")

    for dataset in datasets:
        print(f"  Dataset: {dataset.dataset_id}")
        # Gene symbols and GO annotations must belong to the same species.
        gene_sets_dict = gene_sets_by_namespace[(dataset.species.value, namespace)]

        tec_df: pl.DataFrame = dataset.get_standardized_df(root=standardized_folder)

        rnk_df: pl.DataFrame = get_rnk_from_precalulated_source(tec_df)


        gsea_input: GseaInput = GseaInput(
            rnk=rnk_df,
            gene_sets_dict=gene_sets_dict,
            min_size=7,
            max_size=500,
            permutation_num=10000,
        )

        gsea_result: list[GseaResult] = run_gsea_checkpointed(
            gsea_input,
            checkpoint_file=(output_folder / "checkpoints" / f"{dataset.dataset_id}_{namespace}_std.json"),
            dataset_id=dataset.dataset_id,
            species=dataset.species.value,
        )

        analysed_results.append(
            AnalysedGseaResult(
                gsea_result=gsea_result,
                dataset_id=dataset.dataset_id,
            )
        )

    output_csv: Path = output_folder / f"go_{namespace}_gsea.csv"

    save_gsea_result(output_csv, analysed_results)

    output_excel: Path = save_csv_to_excel(output_csv)

    create_excel_table(output_excel)

    print(f"Saved: {output_csv}")

print("Done")
