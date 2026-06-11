import polars as pl

from data.mouse.cleaned_data.sdk.utils.cleanup_utils import is_unique
from data.rat.sdk.helpers import extract_essential_fields,  convert_ens, convert_rgd
from data.shared_sdk.gene_symbol_cleaning_tools import drop_missing_values, collapse_duplicate_genes_by_mean, \
     cleanup_gene_symbol, remove_unsolved_genes

###### PREPARE DATA ######
rgb_conversion_table_file = r"G:\10. ÜLIÕPILASED\Artur\9. Review\analysis\data\rat\annotation_tables\GENES_RAT_for_rgd.txt"
ens_conversion_table_file = r"G:\10. ÜLIÕPILASED\Artur\9. Review\analysis\data\rat\annotation_tables\mart_export_for_ens.txt"

mtec_lo_file: str = r"mtec_high.csv"
mtec_hi_file: str = r"mtec_low.csv"

current_analyse: str = mtec_hi_file

analyse_file: str = rf"G:\10. ÜLIÕPILASED\Artur\9. Review\analysis\data\rat\raw_data\{current_analyse}"

mtec_df: pl.DataFrame = pl.read_csv(analyse_file)
total_gene_count: int = mtec_df.height

########################
###### DROP NULLS ######
mtec_df, dropped_nulls_count = drop_missing_values(mtec_df)

#################################
###### EXTRACT WHAT NEEDED ######
mtec_df = extract_essential_fields(mtec_df, ["gene_symbol","ctrl1","ctrl2","ctrl3","ko1","ko2","ko3"])

##############################
##### CLEAN GENE SYMBOLS #####
cleaned_symbols: list[str] = []
is_unique_list: list[bool] = []

for gene in mtec_df["gene_symbol"]:
    cleaned_gene: str = cleanup_gene_symbol(gene)

    cleaned_symbols.append(cleaned_gene)
    is_unique_list.append(is_unique(cleaned_gene))

mtec_df = mtec_df.with_columns(
    pl.Series("gene_symbol", cleaned_symbols),
    pl.Series("is_unique", is_unique_list),
)


mouse_ens_count: int = 0
rat_ens_count: int = 0

for gene in mtec_df["gene_symbol"]:
    if gene.startswith("ENSMUSG"):
        mouse_ens_count += 1
    elif gene.startswith("ENSRNOG"):
        rat_ens_count += 1

#########################
###### CONVERT RGD ######
df_rgb_conversion_table = pl.read_csv(
    rgb_conversion_table_file,
    separator="\t",
    comment_prefix="#",
    infer_schema_length=10000,
)

mtec_df, converted_rgd_count = convert_rgd(mtec_df, df_rgb_conversion_table)


del df_rgb_conversion_table

###### CONVERT ENS ######
df_ens_conversion_table = pl.read_csv(
    ens_conversion_table_file,
    separator="\t",
    comment_prefix="#",
    infer_schema_length=10000,
)

mtec_df, converted_ens_count = convert_ens(mtec_df, df_ens_conversion_table)

del df_ens_conversion_table

#################################
###### COLLAPSE DUPLICATES ######
#################################
mtec_df, dropped_duplicates_count = collapse_duplicate_genes_by_mean(mtec_df, 0.1)

mtec_df, unresolved_identifier_count = remove_unsolved_genes(mtec_df)



###### SAVE ######


cleaned_gene_count: int = mtec_df.height

mtec_df.write_csv(rf"G:\10. ÜLIÕPILASED\Artur\9. Review\analysis\data\rat\cleaned_data\cleaned_{current_analyse}")
print(f"Dropped nulls : {dropped_nulls_count}")
print(f"Dropped duplicated genes : {dropped_duplicates_count}")
print(f"Converted RGD IDs: {converted_rgd_count}")
print(f"Converted ENS IDs: {converted_ens_count}")
print(f"Genes after cleaning: {cleaned_gene_count}")

print(f"Total genes: {total_gene_count}")
print(f"Removed unresolved ENS identifiers: {unresolved_identifier_count}")