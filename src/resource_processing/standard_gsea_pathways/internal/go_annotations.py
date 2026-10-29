from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import polars as pl

from src.common.file_io import read_gaf_gz_file
from src.resource_processing.standard_gsea_pathways.internal.models import GoPathway

NOT_QUALIFIER_PATTERN = r"(^|\|)NOT($|\|)"

def prepare_annotations(gaf_data: pl.DataFrame) -> pl.DataFrame:
    required_columns = {
        "go_term",
        "gene_symbol",
        "qualifier",
    }

    missing_columns = required_columns - set(gaf_data.columns)

    if missing_columns:
        raise ValueError(f"GAF table is missing required columns: {sorted(missing_columns)}")

    return (
        gaf_data
        .select(
            pl.col("go_term")
            .cast(pl.String)
            .str.strip_chars(),

            pl.col("gene_symbol")
            .cast(pl.String)
            .str.strip_chars(),

            pl.col("qualifier")
            .cast(pl.String)
            .fill_null(""),
        )
        .filter(
            pl.col("go_term").is_not_null()
            & (pl.col("go_term") != "")
            & pl.col("gene_symbol").is_not_null()
            & (pl.col("gene_symbol") != "")
        )
        .filter(
            ~pl.col("qualifier").str.contains(
                NOT_QUALIFIER_PATTERN
            )
        )
        .select(
            "go_term",
            "gene_symbol",
        )
        .unique()
    )


def find_all_ancestors(go_term: str, pathways_by_id: dict[str, GoPathway], cache: dict[str, frozenset[str]],  active_terms: set[str] | None = None) -> frozenset[str]:
    cached = cache.get(go_term)

    if cached is not None:
        return cached

    if active_terms is None:
        active_terms = set()

    if go_term in active_terms:
        raise ValueError(f"Cycle detected in the GO hierarchy at {go_term}")

    pathway = pathways_by_id.get(go_term)

    if pathway is None:
        return frozenset()

    active_terms.add(go_term)

    ancestors: set[str] = set()

    for parent_id in pathway.parent_go_terms:
        parent = pathways_by_id.get(parent_id)

        if parent is None:
            continue

        # Prevent accidental cross-ontology propagation.
        if parent.namespace != pathway.namespace:
            continue

        ancestors.add(parent_id)

        ancestors.update(
            find_all_ancestors(
                go_term=parent_id,
                pathways_by_id=pathways_by_id,
                cache=cache,
                active_terms=active_terms,
            )
        )

    active_terms.remove(go_term)

    result = frozenset(ancestors)
    cache[go_term] = result

    return result


def attach_genes_to_pathways(pathways: list[GoPathway], annotation_table_file: Path) -> list[GoPathway]:
    annotation_table = read_gaf_gz_file(annotation_table_file)

    annotations = prepare_annotations(annotation_table)

    pathways_by_id: dict[str, GoPathway] = {
        pathway.go_term: pathway
        for pathway in pathways
    }

    ancestor_cache: dict[str, frozenset[str]] = {}

    # Calculate the ancestor closure once for every GO term.
    for go_term in pathways_by_id:
        find_all_ancestors(
            go_term=go_term,
            pathways_by_id=pathways_by_id,
            cache=ancestor_cache,
        )

    genes_by_term: dict[str, set[str]] = {
        go_term: set()
        for go_term in pathways_by_id
    }

    unknown_annotation_count = 0

    # Iterate over unique direct GO-to-gene associations.
    for go_term, gene_symbol in annotations.iter_rows():
        if go_term not in pathways_by_id:
            unknown_annotation_count += 1
            continue

        target_terms = {
            go_term,
            *ancestor_cache[go_term],
        }

        for target_term in target_terms:
            genes_by_term[target_term].add(gene_symbol)

    mapped_pathways = [
        replace(
            pathway,
            genes=sorted(
                genes_by_term[pathway.go_term]
            ),
        )
        for pathway in pathways
    ]

    print(f"Annotations referring to absent or obsolete GO terms: {unknown_annotation_count}")

    return mapped_pathways


def remove_empty_gene_pathways(pathways: list[GoPathway], species: str) -> list[GoPathway]:
    cleaned_pathways = [
        pathway
        for pathway in pathways
        if pathway.genes
    ]

    removed_count = len(pathways) - len(cleaned_pathways)

    print(
        f"{species}: removed {removed_count} pathways without genes; "
        f"kept {len(cleaned_pathways)} pathways with genes."
    )

    return cleaned_pathways


def go_pathways_to_long_dataframe(pathways: list[GoPathway], species: str, include_empty_pathways: bool = False) -> pl.DataFrame:
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

    return pl.DataFrame(
        rows,
        schema={
            "species": pl.String,
            "namespace": pl.String,
            "go_term": pl.String,
            "name": pl.String,
            "gene_symbol": pl.String,
        },
    )