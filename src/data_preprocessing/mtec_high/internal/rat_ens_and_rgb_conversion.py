import polars as pl
from polars import DataFrame

from src.data_preprocessing.mtec_high.internal.gene_symbol_cleaning_tools import is_rgd_id, is_ens_gene_id


def find_ens_gene(conversion_table: pl.DataFrame, query_id: str) -> str | None:
    required_columns: set[str] = {
        "Gene stable ID",
        "Gene name"
    }

    missing_columns: set[str] = required_columns - set(conversion_table.columns)


    if missing_columns:
        raise ValueError(
            "Ensembl conversion table is missing required columns: "
            f"{sorted(missing_columns)}. "
            f"Available columns: {conversion_table.columns}"
        )

    matches: pl.DataFrame = (
        conversion_table
        .filter(
            pl.col("Gene stable ID").cast(pl.String) == query_id
        )
        .select(
            pl.col("Gene name")
            .drop_nulls()
            .unique()
        )
    )

    if matches.height == 0:
        return None

    if matches.height > 2:
        raise ValueError(f"Multiple gene_symbol matches found for {query_id!r}")

    return matches.item(0, "Gene name")

def find_rgd_gene(conversion_table: pl.DataFrame, query_id: str) -> str | None:
    numeric_id: str = query_id.removeprefix("RGD")

    matches: pl.DataFrame = conversion_table.filter(
        pl.col("GENE_RGD_ID").cast(pl.Utf8) == numeric_id
    )

    if matches.height == 0:
        return None

    return matches.item(0, "SYMBOL")

def convert_ens(mtec_df: DataFrame, df_ens_conversion_table: DataFrame) -> tuple[pl.DataFrame, int]:
    converted_ens_count: int = 0
    new_gene_symbols: list[str] = []

    for row in mtec_df.iter_rows(named=True):
        gene_symbol: str = row["gene_symbol"]

        parts: list[str] = [
            part.strip()
            for part in gene_symbol.split("///")
        ]

        for part in parts:
            if is_ens_gene_id(part):
                converted: str | None = find_ens_gene(df_ens_conversion_table, part)

                if converted is not None:
                    gene_symbol = converted
                    converted_ens_count += 1

        new_gene_symbols.append(gene_symbol)

    mtec_df = mtec_df.with_columns(
        pl.Series("gene_symbol", new_gene_symbols)
    )

    return mtec_df, converted_ens_count


def convert_rgd(dataset: DataFrame, df_rgb_conversion_table: DataFrame) -> tuple[pl.DataFrame, int]:
    converted_rgd_count: int = 0
    new_gene_symbols: list[str] = []

    for row in dataset.iter_rows(named=True):
        gene_symbol: str = row["gene_symbol"]

        parts: list[str] = [
            part.strip()
            for part in gene_symbol.split("///")
        ]

        for part in parts:
            if is_rgd_id(part):
                converted: str | None = find_rgd_gene(df_rgb_conversion_table, part)

                if converted is not None:
                    gene_symbol = converted
                    converted_rgd_count += 1

        new_gene_symbols.append(gene_symbol)

    dataset = dataset.with_columns(
        pl.Series("gene_symbol", new_gene_symbols)
    )

    return dataset, converted_rgd_count
