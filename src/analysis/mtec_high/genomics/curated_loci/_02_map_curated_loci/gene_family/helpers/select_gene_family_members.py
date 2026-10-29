from re import Pattern

import polars as pl

def select_gene_family_members(tec_df: pl.DataFrame, gene_patterns: tuple[Pattern[str], ...]) -> pl.DataFrame:
    if not gene_patterns:
        raise ValueError("No gene-family patterns were provided.")

    pattern_strings: list[str] = [
        pattern.pattern
        for pattern in gene_patterns
    ]

    combined_pattern: str = "(?i)(?:" + "|".join(pattern_strings) + ")"

    return tec_df.filter(pl.col("gene_symbol").str.contains(combined_pattern).fill_null(False))
