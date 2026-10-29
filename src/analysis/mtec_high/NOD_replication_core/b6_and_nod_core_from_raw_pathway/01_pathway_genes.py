from pathlib import Path

import polars as pl

from src.analysis.mtec_high.NOD_replication_core.config import CORE_PATHWAYS
from src.common.paths import PROJECT_ROOT

# ------------------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------------------

output_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/NOD_replication_core/b6_and_nod_core_from_raw_pathway/source/01_pathway_genes.csv"

# ------------------------------------------------------------------------------------------
# Extract all mouse genes belonging to the core pathways
# ------------------------------------------------------------------------------------------
pathways_df = pl.read_csv(PROJECT_ROOT / "resources/go/go_pathway_genes.csv").filter(pl.col("species").str.strip_chars().str.to_lowercase() == "mouse")

genes: set[str] = set()

for namespace, pathways in CORE_PATHWAYS.items():
    for pathway in pathways:

        go_term = pathway.split(" | ", 1)[0].strip()

        pathway_df = pathways_df.filter((pl.col("namespace") == namespace.value) & (pl.col("go_term") == go_term))

        if pathway_df.is_empty():
            raise ValueError(f"No mouse pathway genes found for {namespace.value}: {pathway}")

        genes.update(
            gene.strip().lower()
            for gene in pathway_df.get_column("gene_symbol")
            if gene is not None and gene.strip()
        )

# ------------------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------------------

output_df = pl.DataFrame({"leading_gene": sorted(genes)}, schema={"leading_gene": pl.String})
output_file.parent.mkdir(parents=True, exist_ok=True)
output_df.write_csv(output_file)
print(f"Saved {output_df.height} unique mouse pathway genes: {output_file}")
