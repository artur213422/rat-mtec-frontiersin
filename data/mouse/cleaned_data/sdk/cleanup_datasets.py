from pathlib import Path

import polars as pl
from polars import DataFrame

from data.mouse.cleaned_data.sdk.utils.cleanup_utils import is_unique
from data.shared_sdk.gene_symbol_cleaning_tools import drop_missing_values, cleanup_gene_symbol, \
    collapse_duplicate_genes_by_mean, remove_unsolved_genes, cleanup_mouse_rik_genes

print("Running")

folder: str = fr"G:\10. ÜLIÕPILASED\Artur\9. Review\analysis\data\mouse\cleaned_data\converted"
folder_path: Path = Path(folder)

folder: str = fr"G:\10. ÜLIÕPILASED\Artur\9. Review\analysis\data\mouse\cleaned_data\converted"

converted_files: list[Path] = list(folder_path.glob("*_converted.csv"))


for converted_file in converted_files:

    dataset: str = converted_file.stem.removesuffix("_converted")

    raw_data: DataFrame = pl.read_csv(str(converted_file))
    total_row_count: int = raw_data.height

    ## Remove missing values
    raw_data, removed_nulls = drop_missing_values(raw_data)

    ## Clean gene symbols
    cleaned_symbols: list[str] = []
    is_unique_list: list[bool] = []

    for gene in raw_data["gene_symbol"]:
        cleaned_gene: str = cleanup_gene_symbol(gene)

        cleaned_symbols.append(cleaned_gene)
        is_unique_list.append(is_unique(cleaned_gene))

    raw_data = raw_data.with_columns(
        pl.Series("gene_symbol", cleaned_symbols),
        pl.Series("is_unique", is_unique_list),
    )
    ## Remove duplicates
    non_expression_cols: set[str] = {
        "ID_REF",
        "gene_symbol",
        "is_unique",
    }

    raw_data, collapsed_count  = collapse_duplicate_genes_by_mean(raw_data)

    raw_data, unresolved_count  = remove_unsolved_genes(raw_data)

    raw_data, resolved_rik_count = cleanup_mouse_rik_genes(raw_data)

    ## Save
    output_file: str = fr"G:\10. ÜLIÕPILASED\Artur\9. Review\analysis\data\mouse\cleaned_data\cleaned\{dataset}_cleaned.csv"

    unique_gene_count: int = raw_data.height

    raw_data.write_csv(output_file)

    print(f"Dataset: {dataset}")

    print(f"Input rows: {total_row_count}")
    print(f"Dropped null/invalid rows: {removed_nulls}")
    print(f"Removed unresolved genes: {unresolved_count}")
    print(f"Collapsed duplicate rows: {collapsed_count}")
    print(f"Resolved rik genes: {resolved_rik_count}")

    print(f"Unique genes after cleaning: {unique_gene_count}")
    print(f"\n")


