import polars as pl

from src.analysis.mtec_high.simple_clustering.helpers.match_clusters.clean_orthologs import (
    clean_orthologs,
)


def match_rat_clusters(
    matched_mice_df: pl.DataFrame,
    rat_clusters_df: pl.DataFrame,
    reference_dataset: str,
) -> pl.DataFrame:
    """Match complete mouse and rat local clusters by exact ortholog membership."""
    mice = clean_orthologs(matched_mice_df)
    rat = clean_orthologs(rat_clusters_df)

    if rat.filter(pl.col("dataset_id") != "rat_mtec_high").height:
        raise ValueError("Rat input contains another dataset")

    reference = mice.filter(pl.col("dataset_id") == reference_dataset)

    if reference.is_empty():
        raise ValueError(f"Reference dataset '{reference_dataset}' not found")

    def cluster_members(df: pl.DataFrame) -> dict[int, frozenset[str]]:
        clusters = (
            df.filter(pl.col("local_cluster_id").is_not_null())
            .group_by("local_cluster_id")
            .agg(pl.col("ortholog_id").unique().alias("members"))
        )

        return {
            cluster_id: frozenset(members)
            for cluster_id, members in clusters.iter_rows()
            if len(members) > 1
        }

    mouse_clusters = cluster_members(reference)
    rat_clusters = cluster_members(rat)

    mouse_by_members = {
        members: cluster_id
        for cluster_id, members in mouse_clusters.items()
    }

    pairs = [
        (rat_id, mouse_by_members[members])
        for rat_id, members in rat_clusters.items()
        if members in mouse_by_members
    ]

    cluster_map = pl.DataFrame(
        pairs,
        schema={
            "local_cluster_id": rat.schema["local_cluster_id"],
            "matched_cluster_id": reference.schema["local_cluster_id"],
        },
        orient="row",
    )

    # Keep only rat clusters with an exact reference match.
    rat_result = (
        rat.drop("matched_cluster_id", strict=False)
        .join(
            cluster_map,
            on="local_cluster_id",
            how="inner",
            validate="m:1",
        )
    )

    # Use the accepted reference genes to assign all mouse datasets.
    # This also drops mouse genes without an accepted rat counterpart.
    accepted_mouse_ids = cluster_map.get_column("matched_cluster_id")
    reference_map = (
        reference.filter(pl.col("local_cluster_id").is_in(accepted_mouse_ids.to_list()))
        .select(
            "ortholog_id",
            pl.col("local_cluster_id").alias("matched_cluster_id"),
        )
    )

    mice_result = (
        mice.drop("matched_cluster_id", strict=False)
        .join(
            reference_map,
            on="ortholog_id",
            how="inner",
            validate="m:1",
        )
    )

    print(f"Reference mouse clusters (n > 1): {len(mouse_clusters)}")
    print(f"Rat local clusters (n > 1): {len(rat_clusters)}")
    print(f"Strict shared cluster pairs: {len(pairs)}")
    print(f"Unmatched mouse clusters: {len(mouse_clusters) - len(pairs)}")
    print(f"Unmatched rat clusters: {len(rat_clusters) - len(pairs)}")

    return pl.concat([
        mice_result,
        rat_result.select(mice_result.columns),
    ])


def select_shared_orthologs(
    result: pl.DataFrame,
    reference_dataset: str,
) -> pl.DataFrame:
    """Keep genes represented in both members of an accepted cluster pair."""
    reference_keys = (
        result.filter(pl.col("dataset_id") == reference_dataset)
        .select("matched_cluster_id", "ortholog_id")
        .unique()
    )

    rat_keys = (
        result.filter(pl.col("dataset_id") == "rat_mtec_high")
        .select("matched_cluster_id", "ortholog_id")
        .unique()
    )

    shared_keys = reference_keys.join(
        rat_keys,
        on=["matched_cluster_id", "ortholog_id"],
        how="inner",
    )

    return result.join(
        shared_keys,
        on=["matched_cluster_id", "ortholog_id"],
        how="semi",
    )