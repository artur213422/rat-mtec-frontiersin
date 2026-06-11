import polars as pl

from data.mouse.cleaned_data.sdk.utils.cleanup_utils import is_unique
from data.rat.sdk.helpers import extract_essential_fields,  convert_ens, convert_rgd
from data.shared_sdk.gene_symbol_cleaning_tools import drop_missing_values, collapse_duplicate_genes_by_mean, \
    cleanup_gene_symbol, remove_unsolved_genes, cleanup_mouse_rik_genes

###### PREPARE DATA ######

mtec_lo_file: str = r"cleaned_mtec_high.csv"
mtec_hi_file: str = r"cleaned_mtec_low.csv"

current_analyse: str = mtec_hi_file

analyse_file: str = rf"G:\10. ÜLIÕPILASED\Artur\9. Review\analysis\data\rat\cleaned_data\{current_analyse}"

mtec_df: pl.DataFrame = pl.read_csv(analyse_file)
total_gene_count: int = mtec_df.height






mtec_df, converted_rik_count = cleanup_mouse_rik_genes(mtec_df)


mtec_df.write_csv(rf"G:\10. ÜLIÕPILASED\Artur\9. Review\analysis\data\rat\cleaned_data\mousified_{current_analyse}")
