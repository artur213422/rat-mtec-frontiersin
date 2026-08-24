import re
from pathlib import Path
import polars as pl

from src.common.dataset_catalog.api.get_datasets import get_mtec_low
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.feratures.gene_pattern_extraction.methods import find_genes_by_pattern

from src.common.paths import resolve_dir_from_root
from src.common.feratures.gene_pattern_extraction.save_csv import write_report

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


output_dir: Path = resolve_dir_from_root("data/mtec_low/analysis/chemokine_expression_pattern/chemokine_lists")
report_file: Path = output_dir / "chemokine_coverage_report.txt"

datasets: list[DatasetMetadata] = get_mtec_low()

chemokine_patterns: list[re.Pattern[str]] = [
    re.compile(r"^Ccl\d+[a-z]?$", re.IGNORECASE),
    re.compile(r"^Cxcl\d+[a-z]?$", re.IGNORECASE),
    re.compile(r"^Xcl\d+[a-z]?$", re.IGNORECASE),
    re.compile(r"^Cx3cl\d+[a-z]?$", re.IGNORECASE),
]

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

chemokine_dfs: list[pl.DataFrame] = []

for dataset in datasets:
    tec_df: pl.DataFrame = dataset.get_standardized_df()

    tec_df = find_genes_by_pattern(
        data=tec_df,
        gene_patterns=chemokine_patterns,
        gene_col="gene_symbol",
    )

    tec_df = tec_df.with_columns(
        pl.lit(dataset.dataset_id).alias("dataset_id")
    )

    chemokine_dfs.append(tec_df)

chemokines_df: pl.DataFrame = pl.concat(
    chemokine_dfs,
    how="diagonal",
)


##########################################################################################
#                                       SAVE
##########################################################################################

chemokines_df.write_csv(output_dir / "chemokines.csv")

write_report(chemokines_df, report_file, "chemokine")
