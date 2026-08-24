import json
from pathlib import Path

import polars as pl

from src.common.paths import try_find_file_in_project, resolve_dir_from_root
from src.resource_processing.define_tsa_pipeline.internal.tsa_ranking import rank_tsa_group


##########################################################################################
#                                       PREPARE DATA
##########################################################################################

standardized_file: Path = try_find_file_in_project("cleaned_replicate_averaged_tissue_matrix.csv")

tissue_groups_file: Path = try_find_file_in_project("tissue_groups.json")

out_file: Path = resolve_dir_from_root("resources/mouse_gene_atlas/intermediate") / "ranked_tissue_groups.csv"

matrix: pl.DataFrame = pl.read_csv(standardized_file)

tissue_groups: dict[str, list[str]] = json.loads(tissue_groups_file.read_text())

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

ranked_groups: list[pl.DataFrame] = []

for group_name in tissue_groups:
    ranked: pl.DataFrame = rank_tsa_group(
        matrix,
        group_name=group_name,
        tissue_groups=tissue_groups,
    )

    ranked_groups.append(ranked)

ranked_matrix: pl.DataFrame = pl.concat(ranked_groups)

##########################################################################################
#                                           SAVE
##########################################################################################

ranked_matrix.write_csv(out_file)