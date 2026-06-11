from pathlib import Path

import polars as pl
from polars import DataFrame

from data.shared_sdk.gene_symbol_cleaning_tools import is_rgd_id


def find_ens_gene(conversion_table: pl.DataFrame, query_id: str) -> str | None:

    cols_to_search: list[str] = [
        col for col in conversion_table.columns
        if col != "gene_symbol"
    ]

    matches: pl.DataFrame = conversion_table.filter(
        pl.any_horizontal(
            pl.col(col).cast(pl.Utf8) == query_id
            for col in cols_to_search
        )
    )


    if matches.height == 0:
        return None

    if matches.height > 2:
        raise ValueError(f"Multiple gene_symbol matches found for {query_id!r}")


    return matches.item(0, "gene_symbol")



def find_rgd_gene(conversion_table: pl.DataFrame, query_id: str) -> str | None:

    numeric_id: str = query_id.removeprefix("RGD")

    matches: pl.DataFrame = conversion_table.filter(
        pl.col("GENE_RGD_ID").cast(pl.Utf8) == numeric_id
    )

    if matches.height == 0:
        return None

    return matches.item(0, "SYMBOL")



def extract_essential_fields(conversion_table: pl.DataFrame, extract_cols: list[str]) -> pl.DataFrame:
    return conversion_table.select(extract_cols)


def replace_semicolon_gene_splitter(file: str) -> None:
    input_file: Path = Path(file)
    output_file: Path = input_file.with_name(input_file.stem + "_fixed.csv")

    df: DataFrame = pl.read_csv(input_file)

    df = df.with_columns(
        pl.col("gene_symbol")
        .str.replace_all(r"\s*;\s*", " /// ")
        .alias("gene_symbol")
    )

    df.write_csv(output_file)

    print(f"Saved: {output_file}")




def convert_ens(mtec_df: DataFrame, df_ens_conversion_table: DataFrame) -> tuple[pl.DataFrame, int]:
    converted_ens_count: int = 0
    new_gene_symbols: list[str] = []

    for row in mtec_df.iter_rows(named=True):
        gene_symbol: str = row["gene_symbol"]

        parts: list[str] = gene_symbol.split("///")
        for part in parts:
            if is_rgd_id(part):
                converted: str | None = find_ens_gene(df_ens_conversion_table, part)

                if converted is not None:
                    gene_symbol = converted
                    converted_ens_count += 1

        new_gene_symbols.append(gene_symbol)

    mtec_df = mtec_df.with_columns(
        pl.Series("gene_symbol", new_gene_symbols)
    )

    return mtec_df, converted_ens_count


def convert_rgd(mtec_df: DataFrame, df_rgb_conversion_table: DataFrame) -> tuple[pl.DataFrame, int]:
    converted_rgd_count: int = 0
    new_gene_symbols: list[str] = []

    for row in mtec_df.iter_rows(named=True):
        gene_symbol: str = row["gene_symbol"]

        parts: list[str] = gene_symbol.split("///")
        for part in parts:
            if is_rgd_id(part):
                converted: str | None = find_rgd_gene(df_rgb_conversion_table, part)

                if converted is not None:
                    gene_symbol = converted
                    converted_rgd_count += 1

        new_gene_symbols.append(gene_symbol)

    mtec_df = mtec_df.with_columns(
        pl.Series("gene_symbol", new_gene_symbols)
    )

    return mtec_df, converted_rgd_count
