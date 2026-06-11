import polars as pl
from polars import DataFrame

from data.mouse.cleaned_data.sdk.utils.annotations_utils import get_matrix, get_platform_id, get_annotation_table, get_gene_symbol, \
    fix_GSE33878

################# PREPARE DATA #################

series1: str = "GSE85"
series2: str = "GSE2585"
series3: str = "GSE14365"
series4: str = "GSE33878"
on_analyse: str = series4

print(f"Analysing {on_analyse}; that may take a while...")

matrix_file: str = fr"G:\10. ÜLIÕPILASED\Artur\9. Review\analysis\data\mouse\raw_data\{on_analyse}_series_matrix.txt"
output_file: str = fr"G:\10. ÜLIÕPILASED\Artur\9. Review\analysis\data\mouse\cleaned_data\{on_analyse}_converted.csv"
annot_file_folder: str = r"G:\10. ÜLIÕPILASED\Artur\9. Review\analysis\data\mouse\platform_annotation_table"

################# ANALYSE #################
matrix: DataFrame = get_matrix(matrix_file)
platform_id: str = get_platform_id(matrix_file)


annot_table: DataFrame = get_annotation_table(annot_file_folder, platform_id)

gene_symbols: list[str | None] = []


missing_gene_count: int = 0
mapped_gene_count: int = 0

for probe_id in matrix["ID_REF"]:
    gene_name: str | None = get_gene_symbol(str(probe_id), annot_table, "gene_assignment") ## Or Gene Symbol

    if gene_name is None:
        missing_gene_count += 1
        gene_symbols.append(None)
        continue

    mapped_gene_count += 1
    gene_symbols.append(gene_name)

matrix = matrix.with_columns(pl.Series("gene_symbol", gene_symbols))

if on_analyse == series4:
    matrix = fix_GSE33878(matrix)

matrix.write_csv(output_file)
print("Missing gene count: ", missing_gene_count, ", mapped gene count: ", mapped_gene_count)