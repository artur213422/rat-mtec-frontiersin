from pathlib import Path

import polars as pl

from src.common.paths import try_find_file_in_project, resolve_dir_from_root, find_file_in_project
from src.data_preprocessing.mtec_high.internal.orthology_mapping import OrthologRequest, find_ortholog, OrthologSpecies, Status
from src.resource_processing.define_tsa_pipeline.internal.tissue_extraction import extract_list_from_json

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


source_file: Path = try_find_file_in_project("ranked_tissue_groups.csv")
out_file: Path = resolve_dir_from_root("resources/mouse_gene_atlas/intermediate") / "ranked_tissue_groups_with_orthologs.csv"

ortholog_map_file: Path = find_file_in_project("HOM_AllOrganism.rpt.csv")
ortholog_map: pl.DataFrame = pl.read_csv(ortholog_map_file)
ranked_df: pl.DataFrame = pl.read_csv(source_file)


##########################################################################################
#                                       PROCESS DATA
##########################################################################################

genes: list[str] = (
    ranked_df
    .get_column("gene_symbol")
    .unique()
    .to_list()
)

ortholog_ids: list[str | None] = []
has_clean_orthologs: list[bool] = []

for gene in genes:
    result: OrthologRequest = find_ortholog(
        gene,
        OrthologSpecies.MOUSE,
        ortholog_map,
    )

    has_clean_ortholog: bool = True

    match result.status:
        case Status.OK:
            if len(result.rat_symbols) > 1:
                has_clean_ortholog = False

            if len(result.mouse_symbols) > 1:
                has_clean_ortholog = False

        case Status.ORTHOLOGY_ID_NOT_FOUND:
            has_clean_ortholog = False

    ortholog_ids.append(result.ortholog_id)
    has_clean_orthologs.append(has_clean_ortholog)

ortholog_df: pl.DataFrame = pl.DataFrame(
    {
        "gene_symbol": genes,
        "ortholog_id": ortholog_ids,
        "has_clean_ortholog": has_clean_orthologs,
    }
)

ranked_df = ranked_df.join(
    ortholog_df,
    on="gene_symbol",
    how="left",
)

##########################################################################################
#                                           SAVE
##########################################################################################

ranked_df.write_csv(out_file)
