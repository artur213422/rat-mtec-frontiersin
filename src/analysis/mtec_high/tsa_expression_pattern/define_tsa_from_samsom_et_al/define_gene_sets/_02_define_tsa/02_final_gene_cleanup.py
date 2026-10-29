from pathlib import Path
import polars as pl

from src.analysis.mtec_high.tsa_expression_pattern.define_tsa_from_samsom_et_al.define_gene_sets.config.paths import DYNAMIC_STEP_FILE, TSA_FILE
from src.data_preprocessing.mtec_high.internal.gene_symbol_cleaning_tools import drop_invalid_values
from src.analysis.mtec_high.tsa_expression_pattern.define_tsa_from_samsom_et_al.define_gene_sets.config.config import DROP_GENES_WITH_NO_ORTHOLOGS, MAX_SELECTED_TSA_GROUPS
from src.analysis.mtec_high.tsa_expression_pattern.define_tsa_from_samsom_et_al.define_gene_sets.internals.resolve_duplicated_genes import  resolve_duplicated_genes

##########################################################################################
# PREPARE DATA
##########################################################################################

atlas_file: Path = DYNAMIC_STEP_FILE

atlas_df: pl.DataFrame = pl.read_csv(atlas_file)

out_file: Path = TSA_FILE

##########################################################################################
# VALIDATE INPUT
##########################################################################################

required_columns: set[str] = {
    "ID_REF",
    "ortholog_id",
    "gene_symbol",
    "has_clean_ortholog",
    "is_tissue_restricted",
    "primary_j",
    "selected_groups",
}

missing_columns: set[str] = required_columns - set(atlas_df.columns)

if missing_columns:
    raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

##########################################################################################
# CLEAN GENE-SYMBOL ANNOTATION
##########################################################################################

# Clean only the relevant annotation columns.
#
# Do not run drop_invalid_values() on the complete dynamic-step table:
# non-restricted probes legitimately have null boundary/threshold values.
valid_symbol_rows: pl.DataFrame = atlas_df.select("ID_REF", "gene_symbol",)

valid_symbol_rows, dropped_invalid_count = drop_invalid_values(valid_symbol_rows)

print(f"Dropped {dropped_invalid_count} probes with invalid gene-symbol annotations")

# Semi-join retains all dynamic-step columns, but only for probes whose
# gene-symbol annotation passed validation.
atlas_df = atlas_df.join(
    valid_symbol_rows.select("ID_REF"),
    on="ID_REF",
    how="semi",
)

##########################################################################################
# RESTRICT TO CLEAN RAT ORTHOLOGS
##########################################################################################

if DROP_GENES_WITH_NO_ORTHOLOGS:
    before_ortholog_filter: int = atlas_df.height

    atlas_df = atlas_df.filter(
        pl.col("has_clean_ortholog").fill_null(False)
        & pl.col("ortholog_id").is_not_null()
    )

    dropped_ortholog_count: int = (
        before_ortholog_filter - atlas_df.height
    )

    print(f"Dropped {dropped_ortholog_count} probes without a clean rat ortholog")
else:
    # The current resolver groups rows by ortholog_id and therefore
    # cannot process null ortholog IDs.
    if atlas_df["ortholog_id"].null_count() > 0:
        raise ValueError("resolve_duplicated_genes() currently groups by ortholog_id, but the input contains null ortholog IDs")

##########################################################################################
# RESOLVE PROBE-LEVEL CALLS TO GENE-LEVEL CALLS
##########################################################################################

atlas_df = resolve_duplicated_genes(atlas_df=atlas_df, maximum_selected_groups=MAX_SELECTED_TSA_GROUPS)

##########################################################################################
# FINAL OUTPUT COLUMNS
##########################################################################################

if "probe_ids" in atlas_df.columns:
    atlas_df = atlas_df.drop("probe_ids")

##########################################################################################
# SAVE
##########################################################################################

atlas_df.write_csv(out_file)

print(f"Saved: {out_file}")

##########################################################################################
# SUMMARY
##########################################################################################

total_tsa_count: int = atlas_df.height

single_tissue_tsa_count: int = atlas_df.filter(pl.col("primary_j") == 1).height
multi_tissue_tsa_count: int = atlas_df.filter(pl.col("primary_j") > 1).height

print()
print("# FINAL TSA SUMMARY")
print(f"Total gene-level TSAs: {total_tsa_count}")
print(f"Single-tissue TSAs (primary_j = 1): {single_tissue_tsa_count}")
print(f"Multi-tissue TSAs (primary_j > 1): {multi_tissue_tsa_count}")

if total_tsa_count > 0:
    single_tissue_percentage: float = single_tissue_tsa_count / total_tsa_count * 100


    print(f"Single-tissue TSA percentage: {single_tissue_percentage:.2f}%")