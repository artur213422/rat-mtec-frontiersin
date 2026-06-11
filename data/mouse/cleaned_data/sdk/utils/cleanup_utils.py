import polars as pl
from polars import DataFrame

def is_unique(gene_symbol: str) -> bool:
    symbols: list[str] = [
        s.strip()
        for s in gene_symbol.split("///")
        if s.strip() != ""
    ]

    return len(symbols) == 1

