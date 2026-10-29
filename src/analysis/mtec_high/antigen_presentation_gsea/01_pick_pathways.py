import re
from pathlib import Path
import sys

# Support direct execution from any working directory.
PROJECT_ROOT = Path(__file__).resolve().parents[4]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

ANALYSIS_FOLDER = PROJECT_ROOT / "data/mtec_high/derived_data/antigen_presentation_gsea"
STANDARDIZED_FOLDER = PROJECT_ROOT / "data/mtec_high/standardized"

import polars as pl

from src.common.file_io import read_gmt_file
from src.resource_processing.pathways_for_gsea.internal.enums import Namespace

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

output_dir: Path = (ANALYSIS_FOLDER / "pathways")
output_dir.mkdir(parents=True, exist_ok=True)

go_pathways: Path = (PROJECT_ROOT / "resources/go/go_pathway_genes.csv")
reactome_pathways: Path = (PROJECT_ROOT / "resources/reactome/reactome_pathways.gmt")

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

MIN_GENES = 7

##########################################################################################
#                                       FILTER GO PATHWAYS
##########################################################################################

go_pathways_df: pl.DataFrame = pl.read_csv(go_pathways)

go_apc_pathway_genes_df: pl.DataFrame = (
    go_pathways_df
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



go_apc_pathway_genes_df = (
    go_apc_pathway_genes_df
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
#                               FILTER REACTOME PATHWAYS
##########################################################################################
reactome_df: pl.DataFrame = read_gmt_file(reactome_pathways)

pathway_column: str = reactome_df.columns[0]

reactome_apc_pathway_genes_df: pl.DataFrame = (
    reactome_df
    .filter(
        pl.col(pathway_column)
        .fill_null("")
        .str.contains(combined_pattern)
    )
)


##########################################################################################
#                                       SAVE
##########################################################################################

go_output_file: Path = output_dir / "go_antigen_presentation_pathway_genes.csv"
reactome_output_file: Path = output_dir / "reactome_antigen_presentation_pathway_genes.csv"

go_apc_pathway_genes_df.write_csv(go_output_file)
reactome_apc_pathway_genes_df.write_csv(reactome_output_file)

print(f"Saved: {go_output_file}")
print(f"Saved: {reactome_output_file}")

