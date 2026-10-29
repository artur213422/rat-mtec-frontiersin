from pathlib import Path
import sys

# Support direct execution from any working directory.
PROJECT_ROOT = Path(__file__).resolve().parents[7]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.analysis.mtec_high.tsa_expression_pattern.define_tsa_from_samsom_et_al.define_gene_sets.config.paths import (
    DYNAMIC_STEP_FILE,
    NON_TSA_PROBE_FILE,
)

import polars as pl

from src.analysis.mtec_high.tsa_expression_pattern.define_tsa_from_samsom_et_al.define_gene_sets.internals.classify_housekeeping_probes import classify_non_tsa_probes


##########################################################################################
# PREPARE DATA
##########################################################################################

dynamic_step_file: Path = DYNAMIC_STEP_FILE

dynamic_step_df: pl.DataFrame = pl.read_csv(dynamic_step_file)

out_file: Path = NON_TSA_PROBE_FILE
out_file.parent.mkdir(parents=True, exist_ok=True)


##########################################################################################
# PROCESS DATA
##########################################################################################

non_tsa_probe_calls_df: pl.DataFrame = classify_non_tsa_probes(
    dynamic_step_df=dynamic_step_df,
    probe_id_column="ID_REF",
    metadata_columns=[
        "ortholog_id",
        "gene_symbol",
        "has_clean_ortholog",
    ],
)

##########################################################################################
# SUMMARY
##########################################################################################

total_probe_count: int = non_tsa_probe_calls_df.height

tsa_probe_count: int = non_tsa_probe_calls_df.filter(pl.col("is_tsa")).height

non_tsa_probe_count: int = non_tsa_probe_calls_df.filter(pl.col("is_non_tsa")).height

if tsa_probe_count + non_tsa_probe_count != total_probe_count:
    raise ValueError("TSA and non-TSA probe counts do not cover all probes")

print(f"Total probes: {total_probe_count}")
print(f"TSA probes: {tsa_probe_count}")
print(f"Non-TSA probes: {non_tsa_probe_count}")

if total_probe_count > 0:
    print("TSA probe percentage: {tsa_probe_count / total_probe_count * 100:.2f}%")

    print("Non-TSA probe percentage: {non_tsa_probe_count / total_probe_count * 100:.2f}%")

##########################################################################################
# SAVE
##########################################################################################

non_tsa_probe_calls_df.write_csv(out_file)

print(f"Saved: {out_file}")