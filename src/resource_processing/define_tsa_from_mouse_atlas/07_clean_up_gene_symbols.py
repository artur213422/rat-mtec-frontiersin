from pathlib import Path

import polars as pl

from src.common.paths import find_file_in_project
from src.data_preprocessing.mtec_high.internal.gene_symbol_cleaning_tools import unify_gene_symbol_splitters, cleanup_gene_symbol, cleanup_mouse_rik_genes

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

standardized_file: Path = find_file_in_project("gse10246_atlas.csv")
matrix: pl.DataFrame = pl.read_csv(standardized_file)

rik_annotation_table_file: Path = find_file_in_project("MRK_List1.rpt")

##########################################################################################
#                                       PROCESS FILE
##########################################################################################

##  Preparation


### Unify
matrix = unify_gene_symbol_splitters(matrix)

## Clean symbols
cleaned_symbols: list[str] = []

for gene in matrix["gene_symbol"]:
    cleaned_gene: str = cleanup_gene_symbol(gene)
    cleaned_symbols.append(cleaned_gene)

matrix = matrix.with_columns(pl.Series("gene_symbol", cleaned_symbols),)

matrix, resolved_rik_count = cleanup_mouse_rik_genes(matrix, rik_annotation_table_file)


##########################################################################################
#                                       SAV
###########################################################################################

matrix.write_csv(standardized_file)


