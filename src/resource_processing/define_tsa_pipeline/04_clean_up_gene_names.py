from pathlib import Path

import polars as pl

from src.common.paths import resolve_dir_from_root, find_file_in_project, try_find_file_in_project
from src.data_preprocessing.mtec_high.internal.gene_symbol_cleaning_tools import unify_gene_symbol_splitters, drop_invalid_values, \
    cleanup_gene_symbol, cleanup_mouse_rik_genes, generic_collapse_duplicate_genes_by_mean, remove_unsolved_genes
from src.data_preprocessing.mtec_high.internal.models.preprocessing_result.model import PreprocessingResult
from src.data_preprocessing.mtec_high.internal.models.preprocessing_result.ops import save_preprocessing_result

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

standardized_file = try_find_file_in_project("replicate_averaged_tissue_matrix.csv")
matrix: pl.DataFrame = pl.read_csv(standardized_file)

output_folder: Path = resolve_dir_from_root("resources/mouse_gene_atlas/intermediate")
output_file: Path = output_folder /  "cleaned_replicate_averaged_tissue_matrix.csv"

rik_annotation_table_file: Path = find_file_in_project("MRK_List1.rpt")

##########################################################################################
#                                       PROCESS FILE
##########################################################################################

##  Preparation
result: PreprocessingResult = PreprocessingResult("TSA extraction")


result.total_gene_count = matrix.height

### Unify
matrix = unify_gene_symbol_splitters(matrix)


### Drop nulls
matrix, dropped_invalid_count = drop_invalid_values(matrix)

result.dropped_invalid_count = dropped_invalid_count

## Clean symbols
cleaned_symbols: list[str] = []

for gene in matrix["gene_symbol"]:
    cleaned_gene: str = cleanup_gene_symbol(gene)
    cleaned_symbols.append(cleaned_gene)

matrix = matrix.with_columns(
    pl.Series("gene_symbol", cleaned_symbols),
)


matrix, resolved_rik_count = cleanup_mouse_rik_genes(matrix, rik_annotation_table_file)

### Collapse duplicates
matrix, dropped_duplicates_count = generic_collapse_duplicate_genes_by_mean(matrix)
result.dropped_duplicates_count = dropped_duplicates_count

matrix, unresolved_identifier_count = remove_unsolved_genes(matrix)
result.unresolved_identifiers = unresolved_identifier_count

result.cleaned_genes_count = matrix.height


##########################################################################################
#                                       SAV
###########################################################################################

matrix.write_csv(output_file)

save_preprocessing_result(
    result=result,
    output_path=output_folder,
    filename="gene_symbol_cleaning_result.csv",
)


