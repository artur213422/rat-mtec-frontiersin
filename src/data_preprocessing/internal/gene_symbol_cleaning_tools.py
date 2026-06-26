from __future__ import annotations

import re
from pathlib import Path

import polars as pl
from polars import DataFrame

def is_unique(gene_symbol: str) -> bool:
    symbols: list[str] = [
        s.strip()
        for s in gene_symbol.split("///")
        if s.strip() != ""
    ]

    return len(symbols) == 1



def drop_invalid_values(df: DataFrame) -> tuple[DataFrame, int]:
    before: int = df.height

    cleaned = (
        df
        .drop_nulls()
        .filter(
            (pl.col("gene_symbol") != "")
            & (pl.col("gene_symbol") != "---")
            & (~pl.col("gene_symbol").str.to_lowercase().is_in(["abparts", "rnasek"]))
        )
    )

    removed: int = before - cleaned.height

    return cleaned, removed


def unify_gene_symbol_splitters(df: DataFrame) -> DataFrame:
    return df.with_columns(
        pl.col("gene_symbol")
        .str.replace_all(r";|/{2,}", " /// ")
        .alias("gene_symbol")
    )

def is_rik_gene(gene: str) -> bool:
    return bool(re.fullmatch(r"[A-Za-z0-9]+Rik", gene))

def is_rgd_id(gene: str) -> bool:
    return bool(re.fullmatch(r"RGD\d+", gene))


def is_ens_gene_id(gene: str) -> bool:
    return bool(re.fullmatch(r"ENS[A-Z]*G\d+(?:\.\d+)?", gene))

def is_modeled_gene(gene: str) -> bool:
    return bool(re.fullmatch(r"Gm\d+", gene))


def is_loc_gene(gene: str) -> bool:
    return bool(re.fullmatch(r"LOC\d+", gene))


def is_random_text(gene: str) -> bool:
    return gene.lower() == "abparts" or gene.lower() == "rnasek"

def is_identifier_gene(gene: str) -> bool:
    return is_rgd_id(gene) or is_ens_gene_id(gene)


def is_real_gene(gene: str) -> bool:
    return (
        not is_modeled_gene(gene)
        and not is_rik_gene(gene)
        and not is_loc_gene(gene)
        and not is_identifier_gene(gene)
        and not is_random_text(gene)
        and gene != "---"
        and gene != ""
    )

def cleanup_gene_symbol(gene_symbol: str) -> str:
    symbols: list[str] = [
        s.strip()
        for s in gene_symbol.split("///")
        if s.strip() != ""
    ]

    if len(symbols) == 0:
        raise ValueError("Empty symbol list")

    if len(symbols) == 1:
        return symbols[0]

    real_genes: list[str] = [
        gene for gene in symbols
        if is_real_gene(gene)
    ]

    if len(real_genes) == 1:
        return real_genes[0]

    modeled_genes: list[str] = [
        gene for gene in symbols
        if is_modeled_gene(gene)
    ]

    if len(modeled_genes) == 1:
        return modeled_genes[0]

    rik_genes: list[str] = [
        gene for gene in symbols
        if is_rik_gene(gene)
    ]

    if len(rik_genes) == 1:
        return rik_genes[0]

    loc_genes: list[str] = [
        gene for gene in symbols
        if is_loc_gene(gene)
    ]

    identifier_genes: list[str] = [
        gene for gene in symbols
        if is_identifier_gene(gene)
    ]

    if len(identifier_genes) == 1:
        return identifier_genes[0]

    if len(loc_genes) == 1:
        return loc_genes[0]

    abparts_genes: list[str] = [
        gene for gene in symbols
        if is_random_text(gene)
    ]

    if len(abparts_genes) == 1:
        return abparts_genes[0]

    return " /// ".join(symbols)




def collapse_duplicate_genes_by_mean(df: DataFrame, se_floor: float = 0.1,) -> tuple[DataFrame, int]:
    before: int = df.height

    control_cols: list[str] = [
        col
        for col in df.columns
        if col.lower().startswith("ctrl")
    ]

    ko_cols: list[str] = [
        col
        for col in df.columns
        if col.lower().startswith("ko")
    ]

    if len(control_cols) == 0:
        raise ValueError(
            "No control columns found. Expected names like ctrl1, ctrl2, ..."
        )

    if len(ko_cols) == 0:
        raise ValueError(
            "No KO columns found. Expected names like ko1, ko2, ..."
        )

    expression_cols: list[str] = control_cols + ko_cols

    n_control: int = len(control_cols)
    n_ko: int = len(ko_cols)

    df_with_effects: DataFrame = (
        df
        .with_columns(
            [
                pl.col(col).cast(pl.Float64, strict=False).alias(col)
                for col in expression_cols
            ]
        )
        .with_columns(
            pl.mean_horizontal(*control_cols).alias("_control_mean"),
            pl.mean_horizontal(*ko_cols).alias("_ko_mean"),

            pl.concat_list([pl.col(col) for col in control_cols])
            .list.std(ddof=1)
            .alias("_control_sd"),

            pl.concat_list([pl.col(col) for col in ko_cols])
            .list.std(ddof=1)
            .alias("_ko_sd"),
        )
        .with_columns(
            (pl.col("_ko_mean") - pl.col("_control_mean"))
            .alias("_probe_log2_effect"),

            (
                ((pl.col("_control_sd") ** 2) / n_control)
                + ((pl.col("_ko_sd") ** 2) / n_ko)
            )
            .sqrt()
            .alias("_effect_se_raw"),
        )
        .with_columns(
            pl.when(
                pl.col("_effect_se_raw").is_null()
                | pl.col("_effect_se_raw").is_nan()
                | (pl.col("_effect_se_raw") < se_floor)
            )
            .then(pl.lit(se_floor))
            .otherwise(pl.col("_effect_se_raw"))
            .alias("_effect_se")
        )
        .with_columns(
            (1.0 / (pl.col("_effect_se") ** 2)).alias("_effect_weight")
        )
        .with_columns(
            (
                (pl.col("_probe_log2_effect") * pl.col("_effect_weight"))
                .sum()
                .over("gene_symbol")
                /
                pl.col("_effect_weight")
                .sum()
                .over("gene_symbol")
            )
            .alias("_gene_weighted_effect")
        )
        .with_columns(
            (
                (pl.col("_probe_log2_effect") - pl.col("_gene_weighted_effect"))
                / pl.col("_effect_se")
            )
            .alias("_effect_z_residual")
        )
    )

    collapsed: DataFrame = (
        df_with_effects
        .group_by("gene_symbol")
        .agg(
            [pl.len().alias("probe_count")]
            + [
                pl.col(col).mean().alias(col)
                for col in expression_cols
            ]
            + [
                (
                    (pl.col("_effect_z_residual") ** 2)
                    .mean()
                    .sqrt()
                )
                .alias("collapse_inconsistency_z"),

                (
                    pl.col("_probe_log2_effect").max()
                    - pl.col("_probe_log2_effect").min()
                )
                .alias("duplicate_effect_range"),

                pl.col("_probe_log2_effect")
                .mean()
                .alias("mean_probe_log2_effect"),
            ]
        )
        .with_columns(
            pl.when(pl.col("probe_count") > 1)
            .then(pl.col("collapse_inconsistency_z"))
            .otherwise(None)
            .alias("collapse_inconsistency_z"),

            pl.when(pl.col("probe_count") > 1)
            .then(pl.col("duplicate_effect_range"))
            .otherwise(None)
            .alias("duplicate_effect_range"),
        )
    )

    collapsed_genes_count: int = before - collapsed.height

    return collapsed, collapsed_genes_count


def remove_unsolved_genes(data: pl.DataFrame) -> tuple[pl.DataFrame, int]:
    unresolved_identifier_df: pl.DataFrame = data.filter(
        pl.col("gene_symbol").str.starts_with("ENS")
    )

    unresolved_identifier_count: int = unresolved_identifier_df.height


    data = data.filter(
        ~pl.col("gene_symbol").str.starts_with("ENS")
    )

    data = data.filter(
        ~pl.col("gene_symbol").str.contains("///", literal=True)
    )

    data = data.filter(
        ~pl.col("gene_symbol").str.contains(r"^Gm\d+$")
    )

    return data, unresolved_identifier_count


def find_match(df: DataFrame, col_name: str, query: str) -> str | None:
    matches: DataFrame = df.filter(
        pl.col(col_name).cast(pl.Utf8) == query
    )

    if matches.height == 0:
        return None

    return matches.item(0, col_name)



def cleanup_mouse_rik_genes(data: pl.DataFrame, rik_annotation_table_file: Path) -> tuple[pl.DataFrame, int]:
    annotation_table: DataFrame = pl.read_csv(rik_annotation_table_file, separator="\t")

    converted_count: int = 0
    new_gene_symbols: list[str] = []

    for row in data.iter_rows(named=True):
        gene_symbol: str = row["gene_symbol"]

        if is_rik_gene(gene_symbol):
            match: str | None = find_match(
                annotation_table,
                "Marker Symbol",
                gene_symbol
            )

            if match is not None:
                gene_symbol = match
                converted_count += 1

        new_gene_symbols.append(gene_symbol)

    data = data.with_columns(
        pl.Series("gene_symbol", new_gene_symbols)
    )

    return data, converted_count
