import polars as pl


def clean_orthologs(df: pl.DataFrame) -> pl.DataFrame:
    """Keep clean, nonempty single IDs. Never expand composite IDs."""
    result = df.with_columns(pl.col("ortholog_id").cast(pl.String).str.strip_chars())
    result = result.filter(
        pl.col("has_clean_ortholog").fill_null(False)
        & pl.col("ortholog_id").is_not_null()
        & (pl.col("ortholog_id") != "")
        & ~pl.col("ortholog_id").str.contains("///", literal=True)
    )
    duplicates = result.group_by("dataset_id", "ortholog_id").len().filter(pl.col("len") > 1)
    if duplicates.height:
        raise ValueError("Clean ortholog IDs must be unique within each dataset")
    return result
