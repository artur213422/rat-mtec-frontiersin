from pathlib import Path

import polars as pl

from src.analysis.mtec_high.NOD_replication_core.b6_and_node_core_from_leading_genes.config import NOD_DATASETS
from src.analysis.mtec_high.NOD_replication_core.config import CORE_PATHWAYS, GSEA_ANALYSIS
from src.common.paths import PROJECT_ROOT

# ------------------------------------------------------------------------------------------
# Prepare Data
# ------------------------------------------------------------------------------------------

output_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/NOD_replication_core/b6_and_nod_core_from_leading_genes/source/01_nod_pathway_core.csv"

# ------------------------------------------------------------------------------------------
# Extract shared NOD leading-edge genes for each pathway
# ------------------------------------------------------------------------------------------

results: list[dict[str, str]] = []

for namespace, pathways in CORE_PATHWAYS.items():
    if not pathways:
        continue

    gsea_df = pl.read_csv(GSEA_ANALYSIS[namespace])

    for term in pathways:
        pathway_df = gsea_df.filter((pl.col("term") == term) & pl.col("dataset_id").is_in(NOD_DATASETS))

        gene_sets: list[set[str]] = []

        for dataset_id in NOD_DATASETS:
            dataset_df = pathway_df.filter(pl.col("dataset_id") == dataset_id)

            if dataset_df.height != 1:
                raise ValueError(f"Expected one row for {dataset_id}, {term}; found {dataset_df.height}.")

            leading_genes: str | None = dataset_df.item(0, "leading_edge_genes")

            gene_sets.append(
                {
                    gene.strip().lower()
                    for gene in (leading_genes or "").split(";")
                    if gene.strip()
                }
            )

        shared_genes = gene_sets[0] & gene_sets[1]

        for gene in sorted(shared_genes):
            results.append(
                {
                    "namespace": namespace.value,
                    "term": term,
                    "leading_gene": gene,
                }
            )

# ------------------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------------------

output_df = pl.DataFrame(
    results,
    schema={
        "namespace": pl.String,
        "term": pl.String,
        "leading_gene": pl.String,
    },
)

output_file.parent.mkdir(parents=True, exist_ok=True)
output_df.write_csv(output_file)

print(f"Saved {output_df.height} pathway–gene pairs: {output_file}")
