from pathlib import Path

import polars as pl

from src.common.feratures.multi_contditional_limma_voom.induced_specific import (
    join_induction_and_specificity,
)
from src.common.paths import (
    find_file_in_project,
    resolve_dir_from_root,
    resolve_file_from_root,
)
from src.data_preprocessing.mtec_high.internal.orthology_mapping import (
    OrthologSpecies,
    find_ortholog_id,
    map_genes,
)


source_file: Path = resolve_file_from_root("resources/interferon/gse161664_lung/intermediate/ifn_specific_induced_isg/limma_result/ifn_specific.csv")
output_dir: Path = resolve_dir_from_root("resources/interferon/gse161664_lung/standardized/candidates")
output_file: Path = output_dir / "gse161664_ifn_specific_induced_list.csv"



# Read the separate limma contrast results.
limma_df = pl.read_csv(source_file)

# Join target-vs-control and target-vs-other-IFNs-mean by gene.
#
# This produces one row per target IFN and human gene, containing:
# - induction_* statistics
# - specific_* statistics
# - target_condition
source_df = join_induction_and_specificity(
    limma_df,
    target_conditions=("ifna", "ifnb", "ifng", "ifnl"),
)

source_df = source_df.rename({
    "gene_symbol": "human_gene_symbol",
})


# Add mouse ortholog information.
ortholog_map = pl.read_csv(
    find_file_in_project("HOM_AllOrganism.rpt.csv")
)
mapping_rows: list[dict[str, object]] = []

for gene in source_df.get_column(
    "human_gene_symbol"
).unique(maintain_order=True):
    ortholog_id = (
        find_ortholog_id(
            gene,
            OrthologSpecies.HUMAN,
            ortholog_map,
        )
        if gene is not None and gene.strip()
        else None
    )

    mouse_symbols = (
        sorted(
            map_genes(
                ortholog_id,
                OrthologSpecies.MOUSE,
                ortholog_map,
            )
        )
        if ortholog_id is not None
        else []
    )

    mapping_rows.append({
        "human_gene_symbol": gene,
        "ortholog_id": ortholog_id,
        "gene_symbol": ";".join(mouse_symbols) or None,
        "n_mouse_orthologs": len(mouse_symbols),
        "has_clear_ortholog": len(mouse_symbols) == 1,
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

result_df = (
    source_df
    .join(
        mapping_df,
        on="human_gene_symbol",
        how="left",
        validate="m:1",
    )
    .filter(
        pl.col("ortholog_id").is_not_null()
    )
)

# Do not explode ortholog IDs here. That would create multiple rows for
# one IFN–human-gene result.
result_df.write_csv(output_file)

print(f"Saved: {output_file}")