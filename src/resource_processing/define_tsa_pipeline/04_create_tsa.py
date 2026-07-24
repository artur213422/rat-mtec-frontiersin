from pathlib import Path

import polars as pl

from src.common.paths import try_find_file_in_project, resolve_dir_from_root
from src.resource_processing.define_tsa_pipeline.internal.tissue_extraction import  \
    extract_list_from_json
from src.resource_processing.define_tsa_pipeline.internal.tsa_ranking import rank_tsa_candidates, \
    drop_poorest_duplicates

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


standardized_file: Path = resolve_dir_from_root("resources/mouse_gene_atlas/standardized") / "standardized_tissues_matrix.csv"
picked_tissue_file = try_find_file_in_project("picked_tissues.json")

out: Path = resolve_dir_from_root("resources/mouse_gene_atlas/standardized/tsa")

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

tissue_list: list[str] = extract_list_from_json(picked_tissue_file)

matrix: pl.DataFrame = pl.read_csv(standardized_file)

for tissue in tissue_list:
    tsa_df: pl.DataFrame = rank_tsa_candidates(
        matrix=matrix,
        target_tissue=tissue,
        top_n=100,
        min_target_log2=6.0,
        min_specificity_log2=2.0,
    )

    tsa_df = drop_poorest_duplicates(tsa_df)

    if tsa_df.height < 5:
        print(f"{tissue} tissue shape is small: {tsa_df.height}")

    ##########################################################################################
    #                                           SAVE
    ##########################################################################################

    tsa_df.write_csv(out / f"{tissue}_tsa_candidates.csv")

