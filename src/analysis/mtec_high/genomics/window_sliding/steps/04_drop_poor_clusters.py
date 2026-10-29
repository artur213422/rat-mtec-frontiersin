from pathlib import Path

import polars as pl

from src.analysis.mtec_high.genomics.window_sliding.columns import Windows
from src.analysis.mtec_high.genomics.window_sliding.helpers.drop_poor_clusters import drop_poor_clusters
from src.common.excel import create_excel_table
from src.common.file_io import save_csv_to_excel

from src.common.paths import PROJECT_ROOT, find_file_in_project

# -----------------------------------------------------------------------------------------------------------
# Prepare Data
# -----------------------------------------------------------------------------------------------------------

clusters_file: Path = find_file_in_project("03_clusters.parquet")
clusters_df: pl.DataFrame = pl.read_parquet(clusters_file)

output_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/genomics/window_sliding/standardized/clusters.csv"

# -----------------------------------------------------------------------------------------------------------
# Mark Meaningful
# -----------------------------------------------------------------------------------------------------------
datasets: list[str] = clusters_df.get_column("dataset_id").unique().to_list()

marked: list[pl.DataFrame] = []

for dataset in datasets:
    tec_df: pl.DataFrame = clusters_df.filter(pl.col(Windows.DATASET_ID) == dataset)

    tec_df = drop_poor_clusters(tec_df)
    marked.append(tec_df)

# -----------------------------------------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------------------------------------

result: pl.DataFrame = pl.concat(marked)

csv_df: pl.DataFrame = (
    result
    .explode(
        Windows.ORTHOLOG_IDS,
        Windows.GENES,
        Windows.T_STATISTICS,
        Windows.LOG2_FOLD_CHANGES,
        Windows.P_VALUES,
    )
    .rename({
        Windows.ORTHOLOG_IDS: "ortholog_id",
        Windows.GENES: "gene_symbol",
        Windows.T_STATISTICS: "t_statistic",
        Windows.LOG2_FOLD_CHANGES: "log2_fold_change",
        Windows.P_VALUES: "p_val",
    })
    .with_columns(
        pl.col("passenger_genes")
        .list.join("|"),

        pl.col("core_genes")
        .list.join("|"),
    )
)

csv_df.write_csv(output_file)

excel_file: Path = save_csv_to_excel(output_file)
create_excel_table(excel_file)
# -----------------------------------------------------------------------------------------------------------
# Report
# -----------------------------------------------------------------------------------------------------------

summary: pl.DataFrame = (
    result
    .group_by(Windows.DATASET_ID,"aire_direction",)
    .agg(pl.len().alias("cluster_count"))
    .sort( Windows.DATASET_ID, "aire_direction",)
)

print(summary)

