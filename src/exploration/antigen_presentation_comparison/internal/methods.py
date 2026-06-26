from pathlib import Path

import polars as pl



def read_pathways_to_dict(pathways_df: pl.DataFrame) -> dict[str, list[str]]:
    grouped = (
        pathways_df
        .with_columns(
            pathway=pl.concat_str(
                [pl.col("go_term"), pl.lit(" "), pl.col("name")]
            )
        )
        .group_by("pathway")
        .agg(
            pl.col("gene_symbol")
        )
        .sort("pathway")
    )

    return dict(zip(grouped["pathway"], grouped["gene_symbol"], strict=True))