from pathlib import Path

import polars as pl

from src.analysis.mtec_high.top_skin_contributor_genes.config import GO_PATHWAYS, HOLLANDER_PATHWAYS, MICHELSON_PATHWAYS, TSA_PATHWAYS
from src.common.paths import PROJECT_ROOT


def _collect_gsea(
    candidates: pl.DataFrame,
    file: Path,
    pathways: list[str],
    source: str,
    gene_column: str = "leading_edge_genes",
) -> pl.DataFrame:
    """Append one row per leading gene, preserving dataset and pathway provenance."""
    collected = (
        pl.read_csv(file)
        .filter(pl.col("term").is_in(pathways))
        .select(
            pl.col("dataset_id"),
            pl.col(gene_column).str.split(";").alias("leading_gene"),
            pl.lit(source).alias("source"),
            pl.col("term").alias("pathway"),
        )
        .explode("leading_gene")
        .with_columns(pl.col("leading_gene").str.strip_chars())
        .filter(pl.col("leading_gene").is_not_null() & (pl.col("leading_gene") != ""))
    )
    if candidates.width == 0:
        return collected
    return pl.concat([candidates.select(collected.columns), collected], how="vertical_relaxed")


def collect_go_gsea(candidates: pl.DataFrame) -> pl.DataFrame:
    file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/standard_gsea/both_sided_gsea/go_biological_process_gsea.csv"
    return _collect_gsea(candidates, file, GO_PATHWAYS, "GO:BP")


def collect_tsa_gsea(candidates: pl.DataFrame) -> pl.DataFrame:
    file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/tsa_expression/samsom_et_al_tsa_analysis/gsea_per_tissue/gsea_per_tissue.csv"
    return _collect_gsea(candidates, file, TSA_PATHWAYS, "TSA")


def collect_hollander_gsea(candidates: pl.DataFrame) -> pl.DataFrame:
    file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/hollander_thymus_cell_types_gsea/cell_type_gsea.csv"
    return _collect_gsea(candidates, file, HOLLANDER_PATHWAYS, "Hollander et al")


def collect_michelson_gsea(candidates: pl.DataFrame) -> pl.DataFrame:
    file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/mimetic_cells/gsea_analysis.csv"
    return _collect_gsea(candidates, file, MICHELSON_PATHWAYS, "michelson et al", "leading_edge_gene_symbols")
