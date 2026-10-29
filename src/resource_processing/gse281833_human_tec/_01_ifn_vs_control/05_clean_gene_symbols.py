from pathlib import Path
import polars as pl
from src.common.paths import resolve_dir_from_root, find_file_in_project
from src.data_preprocessing.mtec_high.internal.orthology_mapping import find_ortholog_id, OrthologSpecies, map_genes

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

source_file: Path = find_file_in_project("gse281833_tec_ifn_vs_all.csv")

out_file: Path = resolve_dir_from_root("resources/interferon/gse28183_tec/intermediate/ifn_vs_all") / "gse281833_tec_cleaned.csv"

ortholog_map = pl.read_csv(find_file_in_project("HOM_AllOrganism.rpt.csv"))

source_df = pl.read_csv(source_file).rename({
    "gene_symbol": "human_gene_symbol",
})

##########################################################################################
# PROCESS DATA
##########################################################################################

mapping_rows = []

# Map each distinct human symbol once.
for gene in source_df["human_gene_symbol"].unique(maintain_order=True):
    ortholog_id = (
        find_ortholog_id(gene, OrthologSpecies.HUMAN, ortholog_map)
        if gene is not None and gene.strip()
        else None
    )

    mouse_symbols = (
        sorted(map_genes(
            ortholog_id,
            OrthologSpecies.MOUSE,
            ortholog_map,
        ))
        if ortholog_id is not None
        else []
    )

    mapping_rows.append({
        "human_gene_symbol": gene,
        "ortholog_id": ortholog_id,
        "gene_symbol": ";".join(mouse_symbols) or None,
        "n_mouse_orthologs": len(mouse_symbols),
        "has_clear_ortholog": len(mouse_symbols) == 1
    })

mapping_df = pl.DataFrame(
    mapping_rows,
    schema={
        "human_gene_symbol": pl.String,
        "ortholog_id": pl.String,
        "gene_symbol": pl.String,
        "n_mouse_orthologs": pl.Int64,
        "has_clear_ortholog": pl.Boolean,
    },
)

result_df = source_df.join(
    mapping_df,
    on="human_gene_symbol",
    how="left",
    validate="m:1",
)

##########################################################################################
# SAVE
##########################################################################################

result_df.write_csv(out_file)
print(f"Saved: {out_file}")
