import re
from pathlib import Path

import polars as pl

from src.common.paths import find_file_in_project, resolve_dir_from_root
from src.reference_processing.internal.enums import Namespace

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

output_dir: Path = resolve_dir_from_root(
    "data/analysis/peptide_hormone_pattern/pathways"
)
output_dir.mkdir(parents=True, exist_ok=True)

go_pathways: Path = find_file_in_project("go_pathway_genes.csv")

search_patterns: list[re.Pattern[str]] = [
    # Core GO peptide/polypeptide hormone terms
    re.compile(r"\bpeptide hormone\b"),
    re.compile(r"\bpolypeptide hormone\b"),

    # Response/signaling terms, e.g. "response to peptide hormone stimulus"
    re.compile(r"\bresponse to peptide hormone stimulus\b"),
    re.compile(r"\bcellular response to peptide hormone stimulus\b"),
    re.compile(r"\bresponse to polypeptide hormone stimulus\b"),
    re.compile(r"\bcellular response to polypeptide hormone stimulus\b"),

    # Production / secretion / regulation
    re.compile(r"\bpeptide hormone secretion\b"),
    re.compile(r"\bregulation of peptide hormone secretion\b"),
    re.compile(r"\bpositive regulation of peptide hormone secretion\b"),
    re.compile(r"\bnegative regulation of peptide hormone secretion\b"),
    re.compile(r"\bpeptide hormone processing\b"),
    re.compile(r"\bpeptide hormone biosynthetic process\b"),
    re.compile(r"\bpeptide hormone metabolic process\b"),

    # Hormone receptor/signaling, still reasonably focused
    re.compile(r"\bhormone-mediated signaling pathway\b"),
    re.compile(r"\bpeptide hormone receptor\b"),
    re.compile(r"\bpeptide hormone receptor binding\b"),

    # Specific peptide-hormone examples useful for thymus/tolerance context
    re.compile(r"\binsulin\b"),
    re.compile(r"\bglucagon\b"),
    re.compile(r"\bsomatostatin\b"),
    re.compile(r"\bgrowth hormone\b"),
    re.compile(r"\bgonadotropin-releasing hormone\b"),
    re.compile(r"\bthyrotropin-releasing hormone\b"),
    re.compile(r"\bprolactin\b"),
    re.compile(r"\boxytocin\b"),
    re.compile(r"\bvasopressin\b"),
    re.compile(r"\bneuropeptide\b"),
]

combined_pattern: str = "(?i)" + "|".join(
    f"(?:{pattern.pattern})" for pattern in search_patterns
)

MIN_GENES = 5

##########################################################################################
#                                       FILTER PATHWAYS
##########################################################################################

pathways_df: pl.DataFrame = pl.read_csv(go_pathways)

picked_pathway_genes_df: pl.DataFrame = (
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


picked_pathway_genes_df = (
    picked_pathway_genes_df
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

output_file: Path = output_dir / "go_peptide_hormone_pathway_genes.csv"

picked_pathway_genes_df.write_csv(output_file)

print(f"Saved: {output_file}")

unique_pathways_df: pl.DataFrame = (
    picked_pathway_genes_df
    .group_by(["species", "go_term", "name"])
    .agg(
        pl.col("gene_symbol").n_unique().alias("n_genes")
    )
    .sort(["species", "go_term"])
)

print("\nPicked peptide hormone pathways:")
for row in unique_pathways_df.iter_rows(named=True):
    print(
        f"{row['species']} | {row['go_term']} | "
        f"{row['name']} | genes={row['n_genes']}"
    )