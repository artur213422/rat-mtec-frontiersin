from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import polars as pl

from src.common.file_io import read_gaf_gz_file
from src.reference_processing.internal.models import GoPathway


def find_go_genes(gaf_data: pl.DataFrame, target_go_term: str) -> list[str]:
    return (
        gaf_data
        .filter(pl.col("go_term") == target_go_term)
        .filter(
            ~pl.col("qualifier")
            .fill_null("")
            .str.contains(r"(^|\|)NOT($|\|)")
        )
        .select("gene_symbol")
        .drop_nulls()
        .unique()
        .sort("gene_symbol")
        .get_column("gene_symbol")
        .to_list()
    )


def attach_genes_to_pathways(
    pathways: list[GoPathway],
    annotation_table_file: Path,
) -> list[GoPathway]:
    annotation_table: pl.DataFrame = read_gaf_gz_file(annotation_table_file)

    mapped_pathways: list[GoPathway] = []

    for pathway in pathways:
        genes: list[str] = find_go_genes(
            gaf_data=annotation_table,
            target_go_term=pathway.go_term,
        )

        mapped_pathways.append(
            replace(pathway, genes=genes)
        )

    return mapped_pathways


def remove_empty_gene_pathways(
    pathways: list[GoPathway],
    species: str,
) -> list[GoPathway]:
    cleaned_pathways: list[GoPathway] = []

    empty_gene_count: int = 0
    non_empty_gene_count: int = 0

    for pathway in pathways:
        if not pathway.genes:
            empty_gene_count += 1
            continue

        cleaned_pathways.append(pathway)
        non_empty_gene_count += 1

    print(
        f"{species}: removed {empty_gene_count} pathways without genes; "
        f"kept {non_empty_gene_count} pathways with genes."
    )

    return cleaned_pathways


def go_pathways_to_long_dataframe(
    pathways: list[GoPathway],
    species: str,
    include_empty_pathways: bool = False,
) -> pl.DataFrame:
    rows: list[dict[str, str | None]] = []

    for pathway in pathways:
        namespace = pathway.namespace.value

        if pathway.genes:
            for gene in pathway.genes:
                rows.append(
                    {
                        "species": species,
                        "namespace": namespace,
                        "go_term": pathway.go_term,
                        "name": pathway.name,
                        "gene_symbol": gene,
                    }
                )

        elif include_empty_pathways:
            rows.append(
                {
                    "species": species,
                    "namespace": namespace,
                    "go_term": pathway.go_term,
                    "name": pathway.name,
                    "gene_symbol": None,
                }
            )

    return pl.DataFrame(rows)

def remove_invalid_pathways(pathways: list[GoPathway], species: str):
    cleaned_pathways: list[GoPathway] = []

    invalid_path_count: int = 0
    valid_path_count: int = 0

    for pathway in pathways:
        if len(pathway.genes) == 0:
            invalid_path_count += 1
            continue
        else:
            cleaned_pathways.append(pathway)
            valid_path_count += 1

    print(f"{species}: removed {invalid_path_count} invalid pathways; valid pathways: {valid_path_count}")
    return cleaned_pathways