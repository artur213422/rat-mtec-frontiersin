from pathlib import Path

import polars as pl

from src.common.paths import try_find_file_in_project, resolve_dir_from_root, find_file_in_project
from src.data_preprocessing.internal.orthology_mapping import OrthologRequest, find_ortholog, OrthologSpecies, Status
from src.resource_processing.define_tsa_pipeline.internal.tissue_extraction import extract_tissue_matrix, \
    extract_list_from_json
from src.resource_processing.define_tsa_pipeline.internal.tsa_ranking import rank_tsa_candidates

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


tsa_folder: Path = resolve_dir_from_root("resources/mouse_gene_atlas/standardized/tsa")


picked_tissues_file: Path = try_find_file_in_project("picked_tissues.json")
tissues: list[str] = extract_list_from_json(picked_tissues_file)

ortholog_map_file: Path = find_file_in_project("HOM_AllOrganism.rpt.csv")
ortholog_map: pl.DataFrame = pl.read_csv(ortholog_map_file)


##########################################################################################
#                                       PROCESS DATA
##########################################################################################

for tissue in tissues:
    file: Path = try_find_file_in_project(f"{tissue}_tsa_candidates.csv")

    tissue_df: pl.DataFrame = pl.read_csv(file)

    ortholog_ids: list[str | None] = []
    has_clean_orthologs: list[bool] = []

    for row in tissue_df.iter_rows(named=True):
        gene: str = row["gene_symbol"]

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

    tissue_df = tissue_df.with_columns(
        pl.Series("ortholog_id", ortholog_ids),
        pl.Series("has_clean_ortholog", has_clean_orthologs),
    )


    ##########################################################################################
    #                                           SAVE
    ##########################################################################################

    tissue_df.write_csv(file)


