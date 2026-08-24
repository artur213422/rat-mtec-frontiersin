from pathlib import Path

import polars as pl

from src.common.paths import try_find_file_in_project, resolve_dir_from_root, find_file_in_project
from src.data_preprocessing.mtec_high.internal.orthology_mapping import OrthologRequest, find_ortholog, OrthologSpecies, Status
from src.resource_processing.define_tsa_pipeline.internal.tissue_extraction import extract_list_from_json

##########################################################################################
#                                       PREPARE DATA
##########################################################################################
source_file: Path = try_find_file_in_project("ranked_tissue_groups_with_orthologs.csv")
out_file: Path = resolve_dir_from_root("resources/mouse_gene_atlas/standardized") / "ranked_tsa.csv"
report_file: Path = resolve_dir_from_root("resources/mouse_gene_atlas/standardized") / "ranked_tsa_report.txt"

ranked_df: pl.DataFrame = pl.read_csv(source_file)

minimum_mean_enrichment_log2: float = 1.0
minimum_max_specificity_log2: float = 0.0

##########################################################################################
#                                       PROCESS DATA
##########################################################################################
counts_before: pl.DataFrame = (
    ranked_df
    .group_by("tissue_group")
    .len()
    .rename({"len": "before"})
)

ranked_df = ranked_df.filter(
    (pl.col("mean_enrichment_log2") >= minimum_mean_enrichment_log2)
    & (pl.col("max_specificity_log2") > minimum_max_specificity_log2)
)

counts_after: pl.DataFrame = (
    ranked_df
    .group_by("tissue_group")
    .len()
    .rename({"len": "after"})
)

report_df: pl.DataFrame = (
    counts_before
    .join(
        counts_after,
        on="tissue_group",
        how="left",
    )
    .with_columns(
        pl.col("after").fill_null(0),
    )
    .with_columns(
        (pl.col("before") - pl.col("after")).alias("removed"),
    )
    .sort("tissue_group")
)

##########################################################################################
#                                           SAVE
##########################################################################################
ranked_df.write_csv(out_file)

report_lines: list[str] = [
    "TSA filtering report",
    "",
    f"Minimum mean_enrichment_log2: {minimum_mean_enrichment_log2}",
    f"Minimum max_specificity_log2: > {minimum_max_specificity_log2}",
    "",
    f"Total candidates before: {counts_before['before'].sum()}",
    f"Total candidates after:  {report_df['after'].sum()}",
    "",
    "Candidates by tissue/group:",
    "",
]

for row in report_df.iter_rows(named=True):
    report_lines.append(
        f"{row['tissue_group']}: "
        f"{row['before']} -> {row['after']} "
    )

report_file.write_text("\n".join(report_lines), encoding="utf-8",)