from pathlib import Path

import polars as pl

from src.common.paths import resolve_dir_from_root

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

intermediate_dir: Path = resolve_dir_from_root("data/mtec_high/analysis/tsa_expression_pattern/tsa_singed_agreement_analysis/tsa_list/intermediate")

source_file: Path = intermediate_dir / "matched_tsa_selected.csv"
cleaned_file: Path = intermediate_dir / "cleaned_matched_tsa.csv"
report_file: Path = intermediate_dir / "cleaned_matched_tsa_report.txt"

tsa_df: pl.DataFrame = pl.read_csv(source_file)

##########################################################################################
#                                       SETTINGS
##########################################################################################

drop_genes_with_no_orthologs: bool = True

##########################################################################################
#                                       REPORT
##########################################################################################

# Rows with no dataset match are counted as "no match".
no_match: pl.Expr = ~pl.col("found_in_dataset").fill_null(False)

# Only matched genes can meaningfully be classified as having/not having
# a clean ortholog.
no_clean_ortholog: pl.Expr = (
    pl.col("found_in_dataset").fill_null(False)
    & ~pl.col("has_clean_ortholog").fill_null(False)
)

if drop_genes_with_no_orthologs:
    kept: pl.Expr = (
        pl.col("found_in_dataset").fill_null(False)
        & pl.col("has_clean_ortholog").fill_null(False)
    )
else:
    kept = pl.col("found_in_dataset").fill_null(False)


report_df: pl.DataFrame = (
    tsa_df
    .group_by("tissue_group")
    .agg(
        pl.len().alias("total_tsas"),

        pl.col("found_in_dataset")
        .fill_null(False)
        .sum()
        .alias("matched_tsas"),

        no_match.sum().alias("dropped_no_match"),

        (
            no_clean_ortholog.sum()
            if drop_genes_with_no_orthologs
            else pl.lit(0)
        ).alias("dropped_no_ortholog"),

        kept.sum().alias("left"),
    )
    .with_columns(
        (
            pl.col("dropped_no_ortholog")
            / pl.col("matched_tsas")
            * 100
        )
        .round(1)
        .alias("dropped_no_ortholog_pct_of_matched")
    )
    .sort("tissue_group")
)
##########################################################################################
#                                       FILTER DATA
##########################################################################################

tsa_df = tsa_df.filter(
    pl.col("found_in_dataset").fill_null(False)
)

if drop_genes_with_no_orthologs:
    tsa_df = tsa_df.filter(
        pl.col("has_clean_ortholog").fill_null(False)
    )

##########################################################################################
#                                           SAVE
##########################################################################################

tsa_df.write_csv(cleaned_file)

with pl.Config(
    tbl_rows=-1,
    tbl_cols=-1,
    tbl_width_chars=-1,
):
    report_file.write_text(
        str(report_df),
        encoding="utf-8",
    )