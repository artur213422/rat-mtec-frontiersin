from pathlib import Path
from typing import Any
import polars as pl

from src.common.io import load_json, read_gaf_gz_file
from src.common.paths import try_find_file_in_project, resolve_dir_from_root
from src.reference_processing.internal.go_annotations import attach_genes_to_pathways, go_pathways_to_long_dataframe, \
    remove_invalid_pathways
from src.reference_processing.internal.parse_go_basic import map_go_pathways, GoPathway

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

go_pathways_file: Path = try_find_file_in_project("go-basic.json")

mouse_gaf: Path = try_find_file_in_project("mgi.gaf.gz")
rat_gaf: Path = try_find_file_in_project("rgd.gaf.gz")

raw_go_pathways: dict[str, Any] = load_json(go_pathways_file)

output_file: Path = resolve_dir_from_root("reference/go") / "go_pathway_genes.csv"
##########################################################################################
#                                       PROCESS DATA
##########################################################################################

go_pathways: list[GoPathway] = map_go_pathways(raw_go_pathways)


print("Mapping mouse")
mouse_go_pathways: list[GoPathway] = attach_genes_to_pathways(go_pathways, mouse_gaf)

mouse_go_pathways = remove_invalid_pathways(mouse_go_pathways, "mouse")


print("Mapping rat")
rat_go_pathways: list[GoPathway] = attach_genes_to_pathways(go_pathways, rat_gaf)

rat_go_pathways = remove_invalid_pathways(rat_go_pathways, "rat")


del go_pathways
##########################################################################################
#                                       SAVE
##########################################################################################

mouse_df = go_pathways_to_long_dataframe(
    pathways=mouse_go_pathways,
    species="mouse",
)

rat_df = go_pathways_to_long_dataframe(
    pathways=rat_go_pathways,
    species="rat",
)

combined_df = pl.concat([mouse_df, rat_df])

combined_df.write_csv(output_file)



