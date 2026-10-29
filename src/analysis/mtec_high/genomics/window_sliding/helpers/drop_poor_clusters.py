import polars as pl

from src.analysis.mtec_high.genomics.window_sliding.columns import Windows
from src.analysis.mtec_high.genomics.window_sliding.config import P_VAL_THRESHOLD, LOG2_FC_THRESHOLD

def drop_poor_clusters(clusters: pl.DataFrame) -> pl.DataFrame:
    return (
        clusters
        .with_columns(
            _deg_count=_deg_count(),
            _annotated_gene_count=_annotated_gene_count(),
        )
        .filter(
            (pl.col("_deg_count") > 0)
            &
            (pl.col("_annotated_gene_count") > 1)
        )
        .drop(
            "_deg_count",
            "_annotated_gene_count",
        )
    )


def _deg_count() -> pl.Expr:
    return (
        pl.struct(
            [
                pl.col(Windows.P_VALUES),
                pl.col(Windows.LOG2_FOLD_CHANGES),
            ]
        )
        .map_elements(
            lambda x: sum(
                (
                    p <= P_VAL_THRESHOLD
                    and abs(fc) >= LOG2_FC_THRESHOLD
                )
                for p, fc in zip(
                    x[Windows.P_VALUES],
                    x[Windows.LOG2_FOLD_CHANGES],
                )
            ),
            return_dtype=pl.Int64,
        )
    )


def _annotated_gene_count() -> pl.Expr:
    return (
        pl.col(Windows.GENES)
        .list.eval(
            (
                ~pl.element()
                .str.contains(
                    r"^(LOC|RGD|Gm|Rik)"
                )
            )
        )
        .list.sum()
    )
