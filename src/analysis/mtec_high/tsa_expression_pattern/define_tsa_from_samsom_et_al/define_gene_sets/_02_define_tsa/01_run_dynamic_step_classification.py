from pathlib import Path

import polars as pl

from src.analysis.mtec_high.tsa_expression_pattern.define_tsa_from_samsom_et_al.define_gene_sets.config.paths import GROUPED_ATLAS_FILE, DYNAMIC_STEP_FILE
from src.analysis.mtec_high.tsa_expression_pattern.define_tsa_from_samsom_et_al.define_gene_sets.config.config import MINIMUM_EXPRESSION_LINEAR, PRIMARY_SELECTION_RULE, MAX_SELECTED_TSA_GROUPS
from src.analysis.mtec_high.tsa_expression_pattern.define_tsa_from_samsom_et_al.define_gene_sets.internals.dynamic_step_classification import perform_dynamic_step

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

atlas_file: Path = GROUPED_ATLAS_FILE
atlas_df: pl.DataFrame = pl.read_csv(atlas_file)

out_file: Path = DYNAMIC_STEP_FILE
out_file.parent.mkdir(parents=True, exist_ok=True)


##########################################################################################
#                                       PROCESS DATA
##########################################################################################

probe_calls_df: pl.DataFrame = perform_dynamic_step(
    atlas_df,
    probe_id_column="ID_REF",
    metadata_columns=[
        "ortholog_id",
        "gene_symbol",
        "has_clean_ortholog",
    ],
    maximum_selected_groups=MAX_SELECTED_TSA_GROUPS,
    minimum_expression=MINIMUM_EXPRESSION_LINEAR,
    selection_rule=PRIMARY_SELECTION_RULE,
    expression_scale="log2",
)


restricted_probe_count: int = probe_calls_df.filter(pl.col("is_tissue_restricted")).height


total_probe_count: int = probe_calls_df.height

print(f"Total probes: {total_probe_count}")
print(f"Tissue-restricted probes: {restricted_probe_count}")
print(f"Restricted probe percentage: {restricted_probe_count / total_probe_count * 100:.2f}%")

##########################################################################################
#                                           SAVE
##########################################################################################

probe_calls_df.write_csv(out_file)