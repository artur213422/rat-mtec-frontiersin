from pathlib import Path

import polars as pl

from src.common.paths import PROJECT_ROOT

# ------------------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------------------

input_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/NOD_replication_core/b6_and_nod_core_from_raw_pathway/source/05_filtered_candidates.csv"
output_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/NOD_replication_core/b6_and_nod_core_from_raw_pathway/standardized/nod_and_b6_core.csv"

p_value_threshold: float = 0.05
fold_change_threshold: float = 0.5

nod_datasets: tuple[str, ...] = (
    "mouse_gse33878_mtec_high",
    "mouse_gse8564_nod_mtec_high",
)

b6_and_rat_datasets: tuple[str, ...] = (
    "mouse_gse85_mtec_high",
    "mouse_gse2585_mtec_high",
    "mouse_gse8564_b6_mtec_high",
    "rat_mtec_high",
)
min_b6_and_rat_support: int = 2

# ------------------------------------------------------------------------------------------
# Require one DEG plus another dataset changing in the same direction in either group
# ------------------------------------------------------------------------------------------

df = pl.read_csv(input_file, schema_overrides={"ortholog_id": pl.String})

is_deg = (
    (pl.col("p_val") <= p_value_threshold)
    & (pl.col("log2_fold_change").abs() >= fold_change_threshold)
).fill_null(False)

def supported_orthologs(datasets: tuple[str, ...], min_support: int) -> pl.DataFrame:
    # Zero and missing effects have no direction. The supporting non-DEG row
    # needs only a matching sign, with no p-value or effect-size threshold.
    directional_rows = df.filter(
        pl.col("dataset_id").is_in(datasets)
        & pl.col("log2_fold_change").is_not_null()
        & (pl.col("log2_fold_change") != 0)
    ).with_columns(
        pl.col("log2_fold_change").sign().alias("_direction"),
        is_deg.alias("_is_deg"),
    )
    # Count distinct datasets so multiple symbols cannot provide extra support.
    support = directional_rows.group_by("ortholog_id", "_direction").agg(
        pl.col("dataset_id").n_unique().alias("dataset_count"),
        pl.col("_is_deg").any().alias("has_deg"),
    )
    return support.filter(
        pl.col("has_deg") & (pl.col("dataset_count") >= min_support)
    ).select("ortholog_id").unique()


retained_ids = pl.concat([
    supported_orthologs(nod_datasets, len(nod_datasets)),
    supported_orthologs(b6_and_rat_datasets, min_b6_and_rat_support),
]).unique()

# BALB/c is a passenger: retain its rows, but never count it toward support.
# Keep non-DEG rows too, so the output retains the full comparison per gene.
core_df = df.join(retained_ids, on="ortholog_id", how="semi").sort(
    "dataset_id", "ortholog_id", "leading_gene"
)

# ------------------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------------------

output_file.parent.mkdir(parents=True, exist_ok=True)
core_df.write_csv(output_file)

print(f"Retained {retained_ids.height} orthologs and {core_df.height} rows")
print(f"Removed {df.get_column('ortholog_id').n_unique() - retained_ids.height} orthologs")
print(f"Saved to {output_file}")
