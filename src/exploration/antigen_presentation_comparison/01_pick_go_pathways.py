import re
from pathlib import Path

import polars as pl

from src.common.paths import find_file_in_project, resolve_dir_from_root
from src.resource_processing.internal.enums import Namespace

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

output_dir: Path = resolve_dir_from_root(
    "data/analysis/antigen_presentation_comparison/pathways"
)
output_dir.mkdir(parents=True, exist_ok=True)

go_pathways: Path = find_file_in_project("go_pathway_genes.csv")

search_patterns: list[re.Pattern[str]] = [
    re.compile(r"\bantigen processing\b", re.IGNORECASE),
    re.compile(r"\bantigen presentation\b", re.IGNORECASE),
    re.compile(r"\bprocessing and presentation\b", re.IGNORECASE),
    re.compile(r"\bpeptide antigen\b", re.IGNORECASE),
    re.compile(r"\bMHC\b", re.IGNORECASE),
    re.compile(r"\bmajor histocompatibility\b", re.IGNORECASE),
    re.compile(r"\bclass I\b.*\bantigen\b", re.IGNORECASE),
    re.compile(r"\bclass II\b.*\bantigen\b", re.IGNORECASE),
    re.compile(r"\bantigen\b.*\bclass I\b", re.IGNORECASE),
    re.compile(r"\bantigen\b.*\bclass II\b", re.IGNORECASE),
]

combined_pattern: str = "(?i)" + "|".join(
    f"(?:{pattern.pattern})" for pattern in search_patterns
)

##########################################################################################
#                                       FILTER PATHWAYS
##########################################################################################

pathways_df: pl.DataFrame = pl.read_csv(go_pathways)

apc_pathway_genes_df: pl.DataFrame = (
    pathways_df
    .filter(pl.col("namespace") == Namespace.BP.value)
    .filter(pl.col("name").fill_null("").str.contains(combined_pattern))
    .select(
        "species",
        "namespace",
        "go_term",
        "name",
        "gene_symbol",
    )
    .drop_nulls("gene_symbol")
    .unique()
    .sort(["species", "go_term", "gene_symbol"])
)

MIN_GENES = 10

apc_pathway_genes_df = (
    apc_pathway_genes_df
    .with_columns(
        pl.col("gene_symbol")
        .n_unique()
        .over(["species", "go_term"])
        .alias("n_genes")
    )
    .filter(pl.col("n_genes") >= MIN_GENES)
    .drop("n_genes")
)

##########################################################################################
#                                       SAVE
##########################################################################################

output_file: Path = output_dir / "go_antigen_presentation_pathway_genes.csv"

apc_pathway_genes_df.write_csv(output_file)

print(f"Saved: {output_file}")

