import polars as pl


def generate_mice_clustering_matching_summary(cluster_result: pl.DataFrame, ortholog_cluster_result: pl.DataFrame) -> str:
    lines = ["Mouse genes assigned individually to reference regions."]
    for label, frame in [("All genes", cluster_result), ("Clean orthologs", ortholog_cluster_result)]:
        lines.append(label)
        for key, df in frame.group_by("dataset_id", maintain_order=True):
            matched = df.filter(pl.col("matched_cluster_id").is_not_null())
            lines.append(f"{key[0]}: {matched.height}/{df.height} genes mapped; {matched.get_column('matched_cluster_id').n_unique()} reference regions")
    return "\n".join(lines) + "\n"
