from pathlib import Path

import polars as pl

from src.common.paths import find_file_in_project, resolve_dir_from_root
from src.exploration.isg_pathways_analyis.gse161664_human_lung_epithelia_cells.internial.map_orthologs import map_orthologs
from src.exploration.isg_pathways_analyis.gse161664_human_lung_epithelia_cells.internial.summarize import summarize

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

ortholog_map_file: Path = find_file_in_project("HOM_AllOrganism.rpt.csv")
ortholog_map: pl.DataFrame = pl.read_csv(ortholog_map_file)

source_file: Path = find_file_in_project("gse281833_mapped_ifn_vs_all.csv")
summary_file: Path = resolve_dir_from_root("data/analysis/isg_pathways_analysis/gse281833_human_tec/ifn_vs_all_genes/isg_pathways") / "summary.txt"

interferon_isg: pl.DataFrame = pl.read_csv(source_file)

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

interferon_isg = map_orthologs(interferon_isg, ortholog_map)
summarize(interferon_isg, summary_file)

##########################################################################################
#                                           SAVE
##########################################################################################

interferon_isg.write_csv(source_file)


