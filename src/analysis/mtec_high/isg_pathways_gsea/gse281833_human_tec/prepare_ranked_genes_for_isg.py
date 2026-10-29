"""Family-level ISG sensitivity analysis; does not change standardized data.

Use BOTH prepare_ranked_genes_for_isg and prepare_isg_gene_sets_for_recovered.
Linked mouse/rat orthology groups are collapsed consistently across species.
This changes the analysis unit from genes to orthology families. Median ranking
statistics are a descriptive scoring rule, not newly fitted moderated t tests.
The same recovery rules apply to the entire ranked background, not only ISGs.
"""

import csv
import math
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import polars as pl

from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.paths import PROJECT_ROOT
from src.analysis.mtec_high.isg_pathways_gsea.shared.gsea_helpers import prepare_isg_gene_sets


# Supported datasets, not a claim of manual validation of every mapping.
curated_datasets: list[str] = [
    "rat_mtec_high", "mouse_gse85_mtec_high", "mouse_gse2585_mtec_high",
    "mouse_gse33878_mtec_high", "mouse_gse8564_b6_mtec_high",
    "mouse_gse8564_nod_mtec_high", "mouse_gse8564_balbc_mtec_high",
]

ORTHOLOGY_FILE = PROJECT_ROOT / "resources/gene_name_conversion/orthologs/HOM_AllOrganism.rpt.csv"


@lru_cache(maxsize=1)
def _mapping():
    """Resolve connected groups without consulting expression or GSEA scores."""
    with ORTHOLOGY_FILE.open(encoding="utf-8-sig", newline="") as handle:
        rows = [r for r in csv.DictReader(handle)
                if r["Common Organism Name"] in {"rat", "mouse, laboratory"}]
    parent = {r["DB Class Key"]: r["DB Class Key"] for r in rows}

    def find(key):
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    by_symbol = defaultdict(set)
    for r in rows:
        by_symbol[r["Common Organism Name"], r["Symbol"].lower()].add(r["DB Class Key"])
    for groups in by_symbol.values():
        ordered = sorted(groups)
        for group in ordered[1:]:
            left, right = find(ordered[0]), find(group)
            parent[max(left, right)] = min(left, right)

    components = defaultdict(set)
    species = defaultdict(set)
    for group in parent:
        components[find(group)].add(group)
    for r in rows:
        species[find(r["DB Class Key"])].add(r["Common Organism Name"])
    group_to_family = {}
    for root, groups in components.items():
        if species[root] != {"rat", "mouse, laboratory"}:
            continue
        family = root if len(groups) == 1 else "family:" + "_".join(sorted(groups))
        group_to_family.update({g: family for g in groups})

    official = {}
    aliases = defaultdict(set)
    for r in rows:
        family = group_to_family.get(r["DB Class Key"])
        if family is None:
            continue
        species_name = r["Common Organism Name"]
        official[species_name, r["Symbol"].lower()] = family
        for alias in (r.get("Synonyms") or "").split("|"):
            if alias.strip():
                aliases[species_name, alias.strip().lower()].add(family)
    return group_to_family, official, dict(aliases)


def prepare_isg_gene_sets_for_recovered(input_file: Path) -> dict[str, list[str]]:
    """Deduplicate pathway members using the SAME family IDs as the ranks.

    Any member's signature membership assigns the whole family to that signature;
    subtype selectivity therefore remains a family-level approximation.
    """
    mapping, _, _ = _mapping()
    return {
        term: sorted({mapping[g] for g in genes if g in mapping})
        for term, genes in prepare_isg_gene_sets(input_file).items()
    }


def audit_isg_recovery(dataset: DatasetMetadata) -> pl.DataFrame:
    """One row per original measurement, including skipped rows and their reason."""
    if dataset.dataset_id not in curated_datasets:
        raise ValueError(f"Dataset not supported by this workflow: {dataset.dataset_id}")
    species_name = {"rat": "rat", "mouse": "mouse, laboratory"}[str(dataset.species)]
    mapping, official, aliases = _mapping()
    source = dataset.get_standardized_df()
    required = {"gene_symbol", "ortholog_id", "has_clean_ortholog", "ranking_statistic"}
    if missing := required - set(source.columns):
        raise ValueError(f"{dataset.dataset_id}: missing columns {sorted(missing)}")
    audit = []
    for row in source.iter_rows(named=True):
        symbol = (row["gene_symbol"] or "").strip()
        original = (row["ortholog_id"] or "").strip()
        score = row["ranking_statistic"]
        family, method = None, "unresolved"
        key = species_name, symbol.lower()
        if symbol and key in official:
            family, method = official[key], "official_symbol"
        elif symbol and len(aliases.get(key, set())) == 1:
            family, method = next(iter(aliases[key])), "unique_resource_alias"
        else:
            ids = {v.strip() for v in original.split("///") if v.strip()}
            families = {mapping[g] for g in ids if g in mapping}
            if ids and ids <= mapping.keys() and len(families) == 1:
                family, method = next(iter(families)), "stored_group"
        # Refuse conflicts rather than silently reassign an existing measurement.
        stored = {mapping[g.strip()] for g in original.split("///") if g.strip() in mapping}
        if family is not None and stored and stored != {family}:
            family, method = None, "conflicting_symbol_and_group"
        finite = score is not None and math.isfinite(float(score))
        retained = bool(symbol and family is not None and finite)
        audit.append({
            "gene_symbol": symbol, "original_ortholog_id": original,
            "ortholog_id": family, "ranking_statistic": float(score) if finite else None,
            "originally_clean": bool(row["has_clean_ortholog"]),
            "mapping_method": method if finite else "invalid_score",
            "retained": retained,
            "recovered": retained and not bool(row["has_clean_ortholog"]),
        })
    return pl.DataFrame(audit, schema_overrides={
        "ortholog_id": pl.String, "ranking_statistic": pl.Float64,
    })


def prepare_ranked_genes_for_isg(dataset: DatasetMetadata) -> tuple[pl.DataFrame, dict[str, str]]:
    """Recover measurable mappings, then give each family exactly one score.

    Exact resource aliases are accepted only when they identify one family.
    No fuzzy matching, score-dependent mapping, imputation, or score duplication.
    Use audit_isg_recovery(dataset) to inspect every mapping decision.
    """
    source = audit_isg_recovery(dataset).filter(pl.col("retained"))
    if source.is_empty():
        raise ValueError(f"{dataset.dataset_id}: no usable mapped measurements")
    if source["gene_symbol"].str.to_lowercase().n_unique() != source.height:
        raise ValueError(f"{dataset.dataset_id}: duplicate gene symbols; resolve before recovery")
    families = source.group_by("ortholog_id").agg(
        pl.col("ranking_statistic").median(),
        pl.col("gene_symbol").sort().str.join("///"),
    ).sort(["ranking_statistic", "ortholog_id"], descending=[True, False])
    ranks = families.select("ortholog_id", "ranking_statistic")
    symbols = dict(families.select("ortholog_id", "gene_symbol").iter_rows())
    print(f"{dataset.dataset_id}: {source.height} measurements -> {ranks.height} families; "
          f"{source['recovered'].sum()} previously excluded measurements recovered")
    return ranks, symbols
