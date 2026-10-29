
import polars as pl

def drop_invalid_chromosomes(tec_df: pl.DataFrame, dataset_name: str) -> pl.DataFrame:
    n_before: int = tec_df.height

    valid_chromosome: pl.Expr = pl.col("chromosome").str.contains(r"^chr(?:[0-9]+|X|Y)$").fill_null(False)

    removed_chromosomes: list[str] = (
        tec_df
        .filter(~valid_chromosome)
        .get_column("chromosome")
        .drop_nulls()
        .unique()
        .sort()
        .to_list()
    )

    tec_df = tec_df.filter(valid_chromosome)

    n_after: int = tec_df.height
    n_dropped: int = n_before - n_after

    dropped_percent: float = (
        (n_dropped / n_before) * 100
        if n_before > 0
        else 0.0
    )

    print(f"{dataset_name}: dropped {n_dropped}/{n_before} genes; ({dropped_percent:.2f}%). ")
    print()

    return tec_df