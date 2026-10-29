
import polars as pl

from src.analysis.mtec_high.top_skin_contributor_genes.config import MIN_MOUSE_CORE_SIZE

def identify_mouse_core(candidates: pl.DataFrame, pathway: str) -> pl.DataFrame:
    if MIN_MOUSE_CORE_SIZE < 1:
        raise ValueError("MIN_MOUSE_DATASETS must be at least 1")

    mouse_candidates = candidates.filter(
        (pl.col("pathway") == pathway)
        & ~pl.col("dataset_id").str.contains("rat", literal=True)
    )
    return mouse_candidates.filter(pl.col("dataset_id").n_unique().over("gene_symbol") >= MIN_MOUSE_CORE_SIZE)
