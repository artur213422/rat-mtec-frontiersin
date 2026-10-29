from pathlib import Path

import polars as pl

from src.common.paths import resolve_dir_from_root, find_file_in_project

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

ifn_vs_all_isg_file: Path = find_file_in_project("gse161664_ifn_vs_control_isgs.csv")
ifn_specific_isg_file: Path = resolve_dir_from_root("resources/interferon/gse161664_lung/standardized/picked_isgs") / "gse161664_ifn_specific_intersection_isgs.csv"

ifn_vs_all_isg_df = pl.read_csv(ifn_vs_all_isg_file)

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

# An IFN-only ISG occurs in exactly one interferon's ISG set.
ifn_only_orthologs_df = (
    ifn_vs_all_isg_df
    .drop_nulls(subset=["ortholog_id"])
    .select("ortholog_id", "target_condition")
    .unique()
    .group_by("ortholog_id")
    .agg(pl.col("target_condition").n_unique().alias("n_interferons"))
    .filter(pl.col("n_interferons") == 1)
    .select("ortholog_id")
)

result_df = (
    ifn_vs_all_isg_df
    .drop_nulls(subset=["ortholog_id"])
    .join(ifn_only_orthologs_df, on="ortholog_id", how="semi")
    .sort(["target_condition", "log2_fold_change"], descending=[False, True])
)

if result_df.is_empty():
    raise ValueError("No IFN-only ISGs were found.")

for row in (
    result_df
    .group_by("target_condition")
    .len(name="n_ifn_only_isgs")
    .sort("target_condition")
    .iter_rows(named=True)
):
    print(f"{row['target_condition']}: {row['n_ifn_only_isgs']} IFN-only ISGs")

result_df.write_csv(ifn_specific_isg_file)
print(f"Saved: {ifn_specific_isg_file}")
