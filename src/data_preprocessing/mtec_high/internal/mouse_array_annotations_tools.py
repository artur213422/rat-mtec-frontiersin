import polars as pl
from polars import DataFrame





def fix_GSE33878(df: DataFrame) -> DataFrame:
    cleaned_symbols: list[str | None] = []

    for gene_symbol in df["gene_symbol"]:
        if gene_symbol is None:
            cleaned_symbols.append(None)
            continue

        parts = [part.strip() for part in gene_symbol.split("//")]

        if len(parts) < 2:
            cleaned_symbols.append(None)
            continue

        cleaned_symbols.append(parts[1])

    return df.with_columns(
        pl.Series("gene_symbol", cleaned_symbols)
    )


def extract_gene_symbols_from_assignment(gene_assignment: str | None) -> str | None:
    if gene_assignment is None or gene_assignment == "":
        return None

    symbols: list[str] = []

    assignments: list[str] = gene_assignment.split("///")

    for assignment in assignments:
        parts: list[str] = [part.strip() for part in assignment.split("//")]

        if len(parts) < 2:
            continue

        gene_symbol: str = parts[1]

        if gene_symbol == "" or gene_symbol == "---":
            continue

        if gene_symbol not in symbols:
            symbols.append(gene_symbol)

    if len(symbols) == 0:
        return None

    return " /// ".join(symbols)


def fix_weird_GPL6246_gene_symbols(data: pl.DataFrame) -> pl.DataFrame:
    return data.with_columns(
        pl.col("gene_symbol")
        .map_elements(
            extract_gene_symbols_from_assignment,
            return_dtype=pl.String,
        )
        .alias("gene_symbol")
    )