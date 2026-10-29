from pathlib import Path

import polars as pl

from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.paths import PROJECT_ROOT

# ----------------------------------------------------------------------------------------------------
# Prepare data
# ----------------------------------------------------------------------------------------------------

output_folder: Path = PROJECT_ROOT / "data/mtec_high/derived_data/NOD_replication_core/b6_and_nod_core_from_leading_genes"

nod_df: pl.DataFrame = pl.read_csv(output_folder / "source/01_nod_pathway_core.csv")
others_df = pl.read_csv(output_folder / "source/02_others_pathway_core.csv")

nod_genes: set[str] = set(nod_df.get_column("leading_gene"))
other_genes: set[str] = set(others_df.get_column("leading_gene"))


# ----------------------------------------------------------------------------------------------------
# Combine core
# ----------------------------------------------------------------------------------------------------
membership = pl.DataFrame([
    {
        "leading_gene": gene,
        "nod_core": gene in nod_genes,
        "b6_rat_core": gene in other_genes,
        "overlap_class": "shared" if gene in nod_genes & other_genes else "nod_only" if gene in nod_genes else "b6_rat_only",
    }
    for gene in sorted(nod_genes | other_genes)
], schema={"leading_gene": pl.String, "nod_core": pl.Boolean, "b6_rat_core": pl.Boolean, "overlap_class": pl.String})

statistic_columns = ["log2_fold_change", "t_statistic", "p_val", "p_adj", "average_log2_expression"]
annotated = []

for dataset in get_mtec_high():
    reference = dataset.get_standardized_df().select(
        pl.col("gene_symbol").str.strip_chars().str.to_lowercase().alias("leading_gene"),
        "gene_symbol", "ortholog_id",
        *[pl.col(column).cast(pl.Float64) for column in statistic_columns],
    )

    # Preserve the original scripts' case-insensitive symbol-based identity.
    # Missing measurements stay null, never zero. Do not filter on ortholog IDs.
    rows = membership.join(reference, on="leading_gene", how="left", validate="1:1").with_columns(pl.lit(dataset.dataset_id).alias("dataset_id"))
    annotated.append(rows)
    print(f"{dataset.dataset_id}: {rows['gene_symbol'].count()}/{membership.height} genes matched")

standardized = output_folder / "standardized"
standardized.mkdir(parents=True, exist_ok=True)
membership.write_csv(output_folder / "source/03_core_membership.csv")

# ----------------------------------------------------------------------------------------------------
# Save
# ----------------------------------------------------------------------------------------------------

pl.concat(annotated).select(
    "dataset_id", "leading_gene", "gene_symbol", "ortholog_id", "nod_core", "b6_rat_core", "overlap_class", *statistic_columns,
).sort("overlap_class", "leading_gene", "dataset_id").write_csv(standardized / "nod_and_b6_core.csv")

for group, count in membership.group_by("overlap_class").len().sort("overlap_class").iter_rows():
    print(f"{group}: {count} genes")
