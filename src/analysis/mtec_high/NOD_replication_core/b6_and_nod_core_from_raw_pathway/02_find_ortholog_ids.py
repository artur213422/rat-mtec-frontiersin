from pathlib import Path

import polars as pl

from src.common.paths import PROJECT_ROOT, find_file_in_project
from src.data_preprocessing.mtec_high.internal.orthology_mapping import OrthologSpecies, find_ortholog_id

# ------------------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------------------
input_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/NOD_replication_core/b6_and_nod_core_from_raw_pathway/source/01_pathway_genes.csv"


output_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/NOD_replication_core/b6_and_nod_core_from_raw_pathway/source/02_pathway_genes_with_ortholog_id.csv"
# ------------------------------------------------------------------------------------------
# Find ortholog IDs
# ------------------------------------------------------------------------------------------

ortholog_map_file: Path = find_file_in_project("HOM_AllOrganism.rpt.csv")
ortholog_map: pl.DataFrame = pl.read_csv(ortholog_map_file, schema_overrides={"DB Class Key": pl.String})

genes_df = pl.read_csv(input_file, schema_overrides={"leading_gene": pl.String})
ortholog_ids: list[str | None] = [
    find_ortholog_id(gene.strip(), OrthologSpecies.MOUSE, ortholog_map)
    if gene is not None and gene.strip() else None
    for gene in genes_df.get_column("leading_gene")
]
genes_df = genes_df.with_columns(pl.Series("ortholog_id", ortholog_ids, dtype=pl.String))

# ------------------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------------------

output_file.parent.mkdir(parents=True, exist_ok=True)
genes_df.write_csv(output_file)
print(f"Saved {genes_df.height} pathway genes ({genes_df['ortholog_id'].count()} with ortholog IDs): {output_file}")

