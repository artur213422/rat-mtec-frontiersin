import polars as pl


def generate_clustering_summary(
    cluster_result: pl.DataFrame,
    ortholog_cluster_result: pl.DataFrame,
    max_gene_gap_kb: int,
) -> str:
    sections: list[str] = [
        "GENOMIC CLUSTERING SUMMARY",
        "=" * 80,
        f"Max gene gap: {max_gene_gap_kb} kb",
        "",
        _generate_dataset_summary("All genes", cluster_result),
        "",
        _generate_dataset_summary("Orthologs", ortholog_cluster_result),
    ]

    return "\n".join(sections)


def _generate_dataset_summary(
    title: str,
    df: pl.DataFrame,
) -> str:
    sections: list[str] = [
        title,
        "=" * len(title),
    ]

    dataset_ids: list[str] = (
        df
        .get_column("dataset_id")
        .unique()
        .sort()
        .to_list()
    )

    for dataset_id in dataset_ids:
        dataset_df: pl.DataFrame = df.filter(
            pl.col("dataset_id") == dataset_id
        )

        sections.extend(
            [
                "",
                _generate_single_dataset_summary(
                    dataset_id,
                    dataset_df,
                ),
            ]
        )

    return "\n".join(sections)


def _generate_single_dataset_summary(
    dataset_id: str,
    df: pl.DataFrame,
) -> str:
    total_genes: int = df.height

    clustered_df: pl.DataFrame = df.filter(
        pl.col("local_cluster_id").is_not_null()
    )

    unclustered_genes: int = (
        df
        .filter(pl.col("local_cluster_id").is_null())
        .height
    )

    cluster_sizes: pl.DataFrame = (
        clustered_df
        .group_by(["chromosome", "local_cluster_id"])
        .agg(
            pl.len().alias("gene_count")
        )
    )

    cluster_count: int = cluster_sizes.height
    clustered_genes: int = clustered_df.height

    clustered_fraction: float = (
        clustered_genes / total_genes
        if total_genes > 0
        else 0.0
    )

    if cluster_count == 0:
        return "\n".join(
            [
                dataset_id,
                "-" * len(dataset_id),
                f"Genes: {total_genes:,}",
                "Clusters: 0",
                f"Unclustered genes: {unclustered_genes:,}",
            ]
        )

    stats: dict[str, float | int] = (
        cluster_sizes
        .select(
            pl.col("gene_count").median().alias("median"),
            pl.col("gene_count").mean().alias("mean"),
            pl.col("gene_count").max().alias("max"),
        )
        .row(0, named=True)
    )

    largest_cluster: dict[str, object] = (
        cluster_sizes
        .sort("gene_count", descending=True)
        .row(0, named=True)
    )

    return "\n".join(
        [
            dataset_id,
            "-" * len(dataset_id),
            f"Genes: {total_genes:,}",
            f"Clusters: {cluster_count:,}",
            f"Genes in clusters: {clustered_genes:,} "
            f"({clustered_fraction:.1%})",
            f"Unclustered genes: {unclustered_genes:,} "
            f"({1.0 - clustered_fraction:.1%})",
            f"Mean cluster size: {float(stats['mean']):.2f}",
            f"Median cluster size: {float(stats['median']):.1f}",
            f"Largest cluster: {int(stats['max'])} genes "
            f"(chr {largest_cluster['chromosome']}, "
            f"cluster {largest_cluster['local_cluster_id']})",
        ]
    )