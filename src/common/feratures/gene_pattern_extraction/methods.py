from typing import Any

import polars as pl
import re



def find_genes_by_pattern(
    data: pl.DataFrame,
    gene_patterns: list[re.Pattern[str]],
    gene_col: str = "gene_symbol",
) -> pl.DataFrame:
    hits: list[dict[str, Any]] = []

    for row in data.iter_rows(named=True):
        raw_gene_symbol: object = row[gene_col]

        if raw_gene_symbol is None:
            continue

        gene_symbol: str = str(raw_gene_symbol).strip()

        if not gene_symbol:
            continue

        is_hit: bool = any(
            pattern.fullmatch(gene_symbol) is not None
            for pattern in gene_patterns
        )

        if not is_hit:
            continue

        hits.append(row)

    if not hits:
        return data.head(0)

    return (
        pl.DataFrame(
            hits,
            schema=data.schema,
        )
        .sort(gene_col)
    )