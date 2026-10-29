


CLUSTER_GAP_KB: int = 100

PHYSICAL_CLUSTER_COLS: list[str] = [
    "gene_symbol",
    "chromosome",
    "strand",
    "start",
    "stop",

    "dataset_id",

    "cluster_id",
    "cluster_gene_count",

    "gap_to_previous_cluster_bp",
    "gap_to_next_cluster_bp",

    "log2_fold_change",
    "t_statistic",
    "p_val",
    "p_adj",

    "ortholog_id",
    "has_clean_ortholog",
]


CURATED_CLUSTER_COLS: list[str] = [
    "locus_name",
    "locus_definition_type",
    "is_locus_member",

    *PHYSICAL_CLUSTER_COLS,
]