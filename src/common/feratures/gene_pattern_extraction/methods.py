import polars as pl
import re

from src.common.feratures.gene_pattern_extraction.models import GeneReference, SharedGeneReferences, \
    GeneReferencePerDataset


def find_gene_references_by_pattern(
    data: pl.DataFrame,
    gene_patterns: list[re.Pattern[str]],
    gene_col: str = "gene_symbol",
    ortholog_col: str = "ortholog_id",
) -> list[GeneReference]:
    hits: list[GeneReference] = []
    seen: set[tuple[str, int]] = set()

    for row in data.iter_rows(named=True):
        raw_gene_symbol: object = row[gene_col]
        raw_ortholog_id: object = row[ortholog_col]

        if raw_gene_symbol is None or raw_ortholog_id is None:
            continue

        gene_symbol: str = str(raw_gene_symbol).strip()

        if gene_symbol == "":
            continue

        is_chemokine: bool = any(
            pattern.fullmatch(gene_symbol) is not None
            for pattern in gene_patterns
        )

        if not is_chemokine:
            continue

        ortholog_id: int = int(raw_ortholog_id)

        key: tuple[str, int] = (gene_symbol, ortholog_id)

        if key in seen:
            continue

        seen.add(key)

        hits.append(
            GeneReference(
                gene_symbol=gene_symbol,
                ortholog_id=ortholog_id,
            )
        )

    return sorted(
        hits,
        key=lambda hit: hit.gene_symbol,
    )



def find_shared_gene_references(mapped_references: list[GeneReferencePerDataset]) -> SharedGeneReferences:
    if not mapped_references:
        raise ValueError("mapped_references is empty. Why?")

    shared_ortholog_ids: set[int] = {
        hit.ortholog_id
        for hit in mapped_references[0].gen_refs
    }

    union_ortholog_ids: set[int] = set(shared_ortholog_ids)

    hit_by_ortholog_id: dict[int, GeneReference] = {}

    for dataset in mapped_references:
        dataset_ortholog_ids: set[int] = {
            hit.ortholog_id
            for hit in dataset.gen_refs
        }

        shared_ortholog_ids.intersection_update(dataset_ortholog_ids)
        union_ortholog_ids.update(dataset_ortholog_ids)

        for hit in dataset.gen_refs:
            hit_by_ortholog_id.setdefault(hit.ortholog_id, hit)

    missing_ortholog_ids: set[int] = union_ortholog_ids - shared_ortholog_ids

    shared_hits: list[GeneReference] = [
        hit_by_ortholog_id[ortholog_id]
        for ortholog_id in sorted(shared_ortholog_ids)
    ]

    missing_hits: list[GeneReference] = [
        hit_by_ortholog_id[ortholog_id]
        for ortholog_id in sorted(missing_ortholog_ids)
    ]

    return SharedGeneReferences(
        missing_refs=missing_hits,
        shared_refs=shared_hits,
    )