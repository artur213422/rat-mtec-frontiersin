from pathlib import Path

import polars as pl



def read_pathways_to_dict(pathways_df: pl.DataFrame) -> dict[str, list[str]]:
    pathway_expression: pl.Expr

    if "go_term" in pathways_df.columns:
        pathway_expression = pl.concat_str(
            [
                pl.col("go_term"),
                pl.lit(" "),
                pl.col("name"),
            ]
        )
    else:
        pathway_expression = pl.col("pathway")

    grouped: pl.DataFrame = (
        pathways_df
        .with_columns(
            pathway=pathway_expression
        )
        .group_by("pathway")
        .agg(
            pl.col("gene_symbol")
        )
        .sort("pathway")
    )

    return dict(
        zip(
            grouped["pathway"],
            grouped["gene_symbol"],
            strict=True,
        )
    )