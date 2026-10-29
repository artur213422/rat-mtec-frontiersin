from pathlib import Path
import re

import polars as pl

from src.analysis.mtec_high.isg_pathways_gsea.shared.settings import expected_induction_threshold, \
    expected_specific_threshold
from src.common.feratures.multi_contditional_limma_voom.run_multi_conditional_limma_voom import (
    MultiConditionalLimmaVoomInput,
    run_multi_conditional_limma_voom,
    write_target_vs_control_and_other_targets_contrasts,
)
from src.common.paths import resolve_dir_from_root, find_file_in_project

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

source_count_file: Path = find_file_in_project("gse161664_counts.csv")
source_count_df: pl.DataFrame = pl.read_csv(source_count_file)

prepared_for_limma_folder: Path = resolve_dir_from_root("resources/interferon/gse161664_lung/intermediate/ifn_specific_induced_isg/prepared_for_limma")

limma_result_folder: Path = resolve_dir_from_root("resources/interferon/gse161664_lung/intermediate/ifn_specific_induced_isg/limma_result")
prepared_for_limma_folder.mkdir(parents=True, exist_ok=True)
limma_result_folder.mkdir(parents=True, exist_ok=True)

SAMPLE_PREFIXES = {
    "ctrl": "ctrl",
    "ifna": "ifna",
    "ifnb": "ifnb",
    "ifng": "ifng",
    "ifnl": "ifnl3",
}

IFN_NAMES = ("ifna", "ifnb", "ifng", "ifnl")


##########################################################################################
#                                       PROCESS DATA
##########################################################################################
def sample_columns(condition: str) -> list[str]:
    prefix = SAMPLE_PREFIXES[condition]
    columns = [
        column
        for column in source_count_df.columns
        if re.fullmatch(rf"{re.escape(prefix)}\d+", column)
    ]
    if len(columns) < 2:
        raise ValueError(f"Expected at least two samples for {condition}; found {columns}")
    return sorted(columns, key=lambda column: int(column[len(prefix):]))


conditions = ("ctrl", *IFN_NAMES)
condition_sample_columns = {
    condition: sample_columns(condition)
    for condition in conditions
}
selected_columns = [
    sample_id
    for condition in conditions
    for sample_id in condition_sample_columns[condition]
]

counts_df = source_count_df.select(
    "gene_symbol",
    *selected_columns,
).unpivot(
    index="gene_symbol",
    on=selected_columns,
    variable_name="sample_id",
    value_name="count",
)

metadata_df = pl.DataFrame({
    "sample_id": selected_columns,
    "experiment": [1] * len(selected_columns),
    "condition": [
        condition
        for condition in conditions
        for _ in condition_sample_columns[condition]
    ],
})

counts_file = prepared_for_limma_folder / "counts.csv"
metadata_file = prepared_for_limma_folder / "metadata.csv"
contrasts_file = prepared_for_limma_folder / "contrasts.csv"
output_file = limma_result_folder / "ifn_specific.csv"

counts_df.write_csv(counts_file)
metadata_df.write_csv(metadata_file)
write_target_vs_control_and_other_targets_contrasts(
    output_file=contrasts_file,
    design_conditions=conditions,
    target_conditions=IFN_NAMES,
)

run_multi_conditional_limma_voom(
    MultiConditionalLimmaVoomInput(
        counts_file=counts_file,
        metadata_file=metadata_file,
        contrasts_file=contrasts_file,
        output_file=output_file,
        induction_log2_fc = expected_induction_threshold,
        specific_log2_fc = expected_specific_threshold
    )
)

print(f"Saved IFN-specific contrasts: {output_file}")
