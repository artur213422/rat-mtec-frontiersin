from __future__ import annotations

from pathlib import Path

import polars as pl

from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.excel import create_excel_table
from src.common.feratures.gsea.checkpoint import run_gsea_checkpointed
from src.common.feratures.gsea.models import AnalysedGseaResult,GseaInput, GseaMode, GseaResult
from src.common.feratures.gsea.rnk import get_rnk_from_precalulated_source
from src.common.feratures.gsea.save_result import save_gsea_result
from src.common.file_io import save_csv_to_excel
from src.common.paths import PROJECT_ROOT


GO_NAMESPACES: tuple[str, ...] = (
    "biological_process",
    "molecular_function",
    "cellular_component",
)

GSEA_MODES: tuple[GseaMode, ...] = (
    GseaMode.Neg,
    GseaMode.Pos,
)

STANDARDIZED_FOLDER: Path = PROJECT_ROOT / "data/mtec_high/standardized"
OUTPUT_FOLDER: Path = PROJECT_ROOT / "data/mtec_high/derived_data/standard_gsea/one_sided_gsea"

GO_PATHWAYS_FILE: Path = PROJECT_ROOT / "resources/go/go_pathway_genes.csv"


def create_go_gene_sets_dict(
    go_df: pl.DataFrame,
    *,
    namespace: str,
    species: str,
) -> dict[str, list[str]]:
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
        go_df.select(
            pl.col("species").cast(pl.String).str.to_lowercase().str.strip_chars(),
            pl.col("namespace").cast(pl.String).str.to_lowercase().str.strip_chars(),
            pl.col("go_term").cast(pl.String).str.strip_chars(),
            pl.col("name").cast(pl.String).str.strip_chars(),
            pl.col("gene_symbol").cast(pl.String).str.strip_chars(),
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
        .unique(subset=["go_term", "gene_symbol"])
    )

    if namespace_df.is_empty():
        raise ValueError(f"No GO pathways found for namespace={namespace!r}, species={species!r}.")

    gene_sets: dict[str, list[str]] = {}

    for row in (
        namespace_df.group_by("go_term", "name")
        .agg(pl.col("gene_symbol").unique().sort().alias("genes"))
        .sort("go_term")
        .iter_rows(named=True)
    ):
        pathway_key: str = f"{row['go_term']} | {row['name']}"
        gene_sets[pathway_key] = [str(gene) for gene in row["genes"]]

    return gene_sets


def save_results(output_csv: Path, results: list[AnalysedGseaResult]) -> None:
    save_gsea_result(output_csv, results)
    output_excel: Path = save_csv_to_excel(output_csv)
    create_excel_table(output_excel)


def main() -> None:
    OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)

    datasets: list[DatasetMetadata] = get_mtec_high()
    go_pathways_df: pl.DataFrame = pl.read_csv(
        GO_PATHWAYS_FILE,
        infer_schema_length=10_000,
    )

    gene_sets_by_species: dict[str, dict[str, dict[str, list[str]]]] = {
        species: {
            f"go_{namespace}_gsea": create_go_gene_sets_dict(
                go_pathways_df,
                namespace=namespace,
                species=species,
            )
            for namespace in GO_NAMESPACES
        }
        for species in {dataset.species.value for dataset in datasets}
    }

    for mode in GSEA_MODES:
        print(f"Running one-sided GSEA mode: {mode.value}")

        results_by_collection: dict[str, list[AnalysedGseaResult]] = {
            f"go_{namespace}_gsea": [] for namespace in GO_NAMESPACES
        }

        for dataset in datasets:
            print(f"  Dataset: {dataset.dataset_id}")
            # No fallback: missing species annotations must stop the run.
            gene_sets_by_collection = gene_sets_by_species[dataset.species.value]

            tec_df: pl.DataFrame = dataset.get_standardized_df(
                root=STANDARDIZED_FOLDER
            )
            rnk_df: pl.DataFrame = get_rnk_from_precalulated_source(tec_df)

            for collection, gene_sets in gene_sets_by_collection.items():
                print(f"    Collection: {collection}")

                gsea_result: list[GseaResult] = run_gsea_checkpointed(
                    GseaInput(
                        rnk=rnk_df,
                        gene_sets_dict=gene_sets,
                        min_size=7,
                        max_size=500,
                        permutation_num=1_000,
                        gsea_mode=mode,
                    ),
                    checkpoint_file=(OUTPUT_FOLDER / "checkpoints" / f"{dataset.dataset_id}_{collection}_{mode.value}.json"),
                    dataset_id=dataset.dataset_id,
                    species=dataset.species.value,
                )

                results_by_collection[collection].append(
                    AnalysedGseaResult(
                        gsea_result=gsea_result,
                        dataset_id=dataset.dataset_id,
                    )
                )

        for collection, results in results_by_collection.items():
            output_csv: Path = OUTPUT_FOLDER / f"{collection}_{mode.value}.csv"
            save_results(output_csv, results)

    print("Done")


if __name__ == "__main__":
    main()
