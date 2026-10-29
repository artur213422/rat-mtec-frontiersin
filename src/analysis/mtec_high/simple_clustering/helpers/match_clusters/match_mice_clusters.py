import polars as pl


def match_mice_clusters(dataset: pl.DataFrame, reference_dataset: str) -> pl.DataFrame:
    """Assign each gene to its reference region; retain original local IDs."""
    reference = dataset.filter(pl.col("dataset_id") == reference_dataset)
    if reference.is_empty():
        raise ValueError(f"Reference dataset '{reference_dataset}' not found")
    if reference.get_column("gene_symbol").n_unique() != reference.height:
        raise ValueError("Reference gene symbols must be unique")
    reference_map = reference.select(
        "gene_symbol", pl.col("local_cluster_id").alias("matched_cluster_id")
    )
    return dataset.join(reference_map, on="gene_symbol", how="left", validate="m:1")
