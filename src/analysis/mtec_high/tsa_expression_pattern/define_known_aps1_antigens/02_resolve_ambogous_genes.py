from pathlib import Path

import polars as pl

from src.common.paths import find_file_in_project, resolve_dir_from_root

##########################################################################################
# PREPARE DATA
##########################################################################################

aps1_antigenes_file: Path = find_file_in_project("aps1_antigens_with_orthologs.csv")

aps1_antigenes_df: pl.DataFrame = pl.read_csv(aps1_antigenes_file)

out_file: Path = (resolve_dir_from_root("data/mtec_high/derived_data/tsa_expression/aps1_known_antigens/source/intermediate/") /
                  "aps1_antigens_with_resolved_orthologs.csv")

##########################################################################################
# RESOLVE
##########################################################################################
# Expand ambiguous one-to-many ortholog mappings into separate rows.
# Example:
#   gene_symbol = "Ins1, Ins2"
#   ortholog_id = "51559129///51559130"

aps1_antigenes_df = aps1_antigenes_df.with_row_index("_row_order")

clear_df = aps1_antigenes_df.filter(pl.col("has_clean_ortholog"))

ambiguous_df = (
    aps1_antigenes_df
    .filter(~pl.col("has_clean_ortholog"))
    .with_columns(
        pl.col("gene_symbol")
        .str.split(",")
        .list.eval(pl.element().str.strip_chars())
        .alias("_gene_symbols"),

        pl.col("ortholog_id")
        .cast(pl.String)
        .str.split("///")
        .list.eval(pl.element().str.strip_chars())
        .alias("_ortholog_ids"),
    )
)

# Validate that each ambiguous row has equally many gene symbols and IDs.
invalid_mappings = ambiguous_df.filter(
    pl.col("_gene_symbols").list.len()
    != pl.col("_ortholog_ids").list.len()
)

if invalid_mappings.height:
    raise ValueError(
        "The number of gene_symbol and ortholog_id values differs after splitting:\n"
        f"{invalid_mappings.select(
            '_row_order',
            'gene_symbol',
            'ortholog_id',
            pl.col('_gene_symbols').list.len().alias('gene_count'),
            pl.col('_ortholog_ids').list.len().alias('ortholog_id_count'),
        )}"
    )

ambiguous_df = (
    ambiguous_df
    # Exploding both columns together pairs values by list position.
    .explode("_gene_symbols", "_ortholog_ids")
    .drop("gene_symbol", "ortholog_id")
    .rename(
        {
            "_gene_symbols": "gene_symbol",
            "_ortholog_ids": "ortholog_id",
        }
    )
)

aps1_antigenes_df = (
    pl.concat(
        [clear_df, ambiguous_df],
        how="diagonal_relaxed",
    )
    .sort("_row_order")
    .drop("_row_order")
)

aps1_antigenes_df = aps1_antigenes_df.drop("has_clean_ortholog")

##########################################################################################
# SAVE
##########################################################################################
aps1_antigenes_df.write_csv(out_file)

