from pathlib import Path
import sys

# Support direct execution from any working directory.
PROJECT_ROOT = Path(__file__).resolve().parents[7]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.analysis.mtec_high.tsa_expression_pattern.define_tsa_from_samsom_et_al.define_gene_sets.config.paths import (
    NON_TSA_PROBE_FILE,
    NON_TSA_FILE,
)

import polars as pl

from src.data_preprocessing.mtec_high.internal.gene_symbol_cleaning_tools import (
    drop_invalid_values,
)
from src.analysis.mtec_high.tsa_expression_pattern.define_tsa_from_samsom_et_al.define_gene_sets.config.config import (
    DROP_GENES_WITH_NO_ORTHOLOGS,
)

##########################################################################################
# PREPARE DATA
##########################################################################################

probe_calls_file: Path = NON_TSA_PROBE_FILE

probe_calls_df: pl.DataFrame = pl.read_csv(
    probe_calls_file
)

out_file: Path = NON_TSA_FILE
out_file.parent.mkdir(parents=True, exist_ok=True)

##########################################################################################
# VALIDATE INPUT
##########################################################################################

required_columns: set[str] = {
    "ID_REF",
    "ortholog_id",
    "gene_symbol",
    "has_clean_ortholog",
    "is_tsa",
    "is_non_tsa",
}

missing_columns: set[str] = (
    required_columns - set(probe_calls_df.columns)
)

if missing_columns:
    raise ValueError(
        f"Missing required columns: {sorted(missing_columns)}"
    )

if probe_calls_df["ID_REF"].null_count() > 0:
    raise ValueError(
        "Input contains null probe IDs"
    )

if probe_calls_df["is_tsa"].null_count() > 0:
    raise ValueError(
        "Input contains null TSA calls"
    )

if probe_calls_df["is_non_tsa"].null_count() > 0:
    raise ValueError(
        "Input contains null non-TSA calls"
    )

duplicated_probe_count: int = probe_calls_df.filter(
    pl.col("ID_REF").is_duplicated()
).height

if duplicated_probe_count > 0:
    raise ValueError(
        f"Input contains {duplicated_probe_count} duplicated probe IDs"
    )

##########################################################################################
# VALIDATE COMPLEMENTARY CALLS
##########################################################################################

invalid_complement_count: int = probe_calls_df.filter(
    pl.col("is_tsa") == pl.col("is_non_tsa")
).height

if invalid_complement_count > 0:
    raise ValueError(
        f"{invalid_complement_count} probes have invalid TSA/non-TSA calls. "
        "Exactly one of is_tsa and is_non_tsa must be True."
    )

##########################################################################################
# CLEAN GENE-SYMBOL ANNOTATION
##########################################################################################

valid_symbol_rows: pl.DataFrame = probe_calls_df.select(
    "ID_REF",
    "gene_symbol",
)

valid_symbol_rows, dropped_invalid_count = (
    drop_invalid_values(
        valid_symbol_rows
    )
)

print(
    f"Dropped {dropped_invalid_count} probes "
    "with invalid gene-symbol annotations"
)

probe_calls_df = probe_calls_df.join(
    valid_symbol_rows.select("ID_REF"),
    on="ID_REF",
    how="semi",
)

##########################################################################################
# RESTRICT TO CLEAN RAT ORTHOLOGS
##########################################################################################

if DROP_GENES_WITH_NO_ORTHOLOGS:
    before_ortholog_filter: int = probe_calls_df.height

    probe_calls_df = probe_calls_df.filter(
        pl.col("has_clean_ortholog").fill_null(False)
        & pl.col("ortholog_id").is_not_null()
    )

    dropped_ortholog_count: int = (
        before_ortholog_filter
        - probe_calls_df.height
    )

    print(
        f"Dropped {dropped_ortholog_count} probes "
        "without a clean rat ortholog"
    )

else:
    if probe_calls_df["ortholog_id"].null_count() > 0:
        raise ValueError(
            "The non-TSA resolver groups probes by ortholog_id, "
            "but the input contains null ortholog IDs"
        )

##########################################################################################
# RESOLVE PROBE CALLS TO GENE CALLS
##########################################################################################

# Gene-level rule:
#
# - all probes TSA:
#       gene is not included in the non-TSA resource;
#
# - all probes non-TSA:
#       gene is included;
#
# - mixed TSA and non-TSA probes:
#       True wins, so the gene is included as non-TSA.

gene_calls_df: pl.DataFrame = (
    probe_calls_df
    .group_by(
        "ortholog_id",
        maintain_order=True,
    )
    .agg(
        pl.col("gene_symbol")
        .sort()
        .first()
        .alias("gene_symbol"),

        pl.col("gene_symbol")
        .n_unique()
        .alias("_gene_symbol_count"),

        pl.len()
        .alias("probe_count"),

        pl.col("is_non_tsa")
        .all()
        .alias("_all_probes_non_tsa"),

        pl.col("is_non_tsa")
        .any()
        .alias("is_non_tsa"),

        pl.col("ID_REF")
        .sort()
        .str.join(";")
        .alias("probe_ids"),
    )
    .with_columns(
        (
            pl.col("_gene_symbol_count") > 1
        ).alias("gene_symbol_conflict"),

        (
            pl.col("_all_probes_non_tsa")
            == pl.col("is_non_tsa")
        ).alias("probe_calls_agree"),
    )
)

##########################################################################################
# REPORT PROBE AGREEMENT
##########################################################################################

total_gene_count: int = gene_calls_df.height

multiple_probe_gene_count: int = gene_calls_df.filter(
    pl.col("probe_count") > 1
).height

unanimous_non_tsa_count: int = gene_calls_df.filter(
    pl.col("_all_probes_non_tsa")
).height

unanimous_tsa_count: int = gene_calls_df.filter(
    ~pl.col("is_non_tsa")
).height

mixed_probe_call_count: int = gene_calls_df.filter(
    pl.col("is_non_tsa")
    & ~pl.col("_all_probes_non_tsa")
).height

print()
print("# PROBE-TO-GENE NON-TSA RESOLUTION")
print()
print(
    f"Total ortholog genes examined: "
    f"{total_gene_count}"
)
print(
    f"Genes represented by multiple probes: "
    f"{multiple_probe_gene_count}"
)
print()
print(
    f"Genes unanimously called non-TSA: "
    f"{unanimous_non_tsa_count}"
)
print(
    f"Genes unanimously called TSA: "
    f"{unanimous_tsa_count}"
)
print(
    f"Genes with mixed TSA/non-TSA probes: "
    f"{mixed_probe_call_count}"
)
print(
    f"Mixed-call genes retained because True wins: "
    f"{mixed_probe_call_count}"
)

##########################################################################################
# RETAIN NON-TSA GENES
##########################################################################################

non_tsa_df: pl.DataFrame = (
    gene_calls_df
    .filter(
        pl.col("is_non_tsa")
    )
    .select(
        "ortholog_id",
        "gene_symbol",
        "is_non_tsa",
        "probe_count",
        "probe_calls_agree",
        "gene_symbol_conflict",
        "probe_ids",
    )
    .sort(
        [
            "probe_calls_agree",
            "gene_symbol",
        ],
        descending=[
            True,
            False,
        ],
    )
)

##########################################################################################
# VALIDATE OUTPUT
##########################################################################################

duplicated_gene_count: int = non_tsa_df.filter(
    pl.col("ortholog_id").is_duplicated()
).height

if duplicated_gene_count > 0:
    raise ValueError(
        f"Output contains {duplicated_gene_count} "
        "duplicated ortholog IDs"
    )

if not non_tsa_df.get_column("is_non_tsa").all():
    raise ValueError(
        "Output contains genes that are not classified as non-TSA"
    )

##########################################################################################
# SAVE
##########################################################################################

non_tsa_df.write_csv(
    out_file
)

print()
print(f"Saved: {out_file}")

##########################################################################################
# FINAL SUMMARY
##########################################################################################

total_non_tsa_count: int = non_tsa_df.height

non_tsa_percentage: float = (
    total_non_tsa_count
    / total_gene_count
    * 100
    if total_gene_count > 0
    else 0.0
)

print()
print("# FINAL NON-TSA SUMMARY")
print(
    f"Final gene-level non-TSAs: "
    f"{total_non_tsa_count}"
)
print(
    f"Non-TSA percentage among examined genes: "
    f"{non_tsa_percentage:.2f}%"
)