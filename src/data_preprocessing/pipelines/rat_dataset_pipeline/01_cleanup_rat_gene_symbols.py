from pathlib import Path
import polars as pl

from src.common.models.dataset_metadata.model import DatasetCatalog, DatasetMetadata
from src.common.models.dataset_metadata.ops import build_dataset_catalog
from src.common.paths import find_file_in_project, resolve_dir_from_root
from src.data_preprocessing.internal.gene_symbol_cleaning_tools import drop_invalid_values, cleanup_gene_symbol, \
    unify_gene_symbol_splitters, collapse_duplicate_genes_by_mean, remove_unsolved_genes, is_unique
from src.data_preprocessing.internal.models.preprocessing_result.model import PreprocessingResult
from src.data_preprocessing.internal.models.preprocessing_result.ops import save_preprocessing_result
from src.data_preprocessing.internal.rat_ens_and_rgb_conversion import convert_rgd, convert_ens

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

ens_conversion_table_file: Path = find_file_in_project("biomart_ens_rat_symbol_annotations.txt")
rgb_conversion_table_file: Path = find_file_in_project("rat_rgb_symbol_annotations.txt")

datasets: DatasetCatalog = build_dataset_catalog()
rat_datasets: list[DatasetMetadata] = datasets.get_rat_datasets()

output: Path = resolve_dir_from_root("data/intermediate/rat")
raw_data_path: Path = resolve_dir_from_root("data/raw/rat")

##########################################################################################
#                                       PROCESS FILE
##########################################################################################
for dataset in rat_datasets:
    ### Preparation
    print("Processing:", dataset.file.name)
    result: PreprocessingResult = PreprocessingResult(dataset.dataset_name)

    source_file: Path = raw_data_path / dataset.file.name

    if not source_file.exists():
        raise FileNotFoundError(f"{source_file} does not exist")

    tec_df: pl.DataFrame = pl.read_csv(source_file)

    result.total_gene_count = tec_df.height

    selected_cols: list[str] = ["gene_symbol","ctrl1","ctrl2","ctrl3","ko1","ko2","ko3"]
    tec_df = tec_df.select(selected_cols)

    ### Unify
    tec_df = unify_gene_symbol_splitters(tec_df)

    ### Drop nulls
    tec_df, dropped_invalid_count = drop_invalid_values(tec_df)

    result.dropped_invalid_count = dropped_invalid_count

    ## Clean symbols
    cleaned_symbols: list[str] = []

    for gene in tec_df["gene_symbol"]:
        cleaned_gene: str = cleanup_gene_symbol(gene)
        cleaned_symbols.append(cleaned_gene)

    tec_df = tec_df.with_columns(
        pl.Series("gene_symbol", cleaned_symbols),
    )

    ## Process RGD ######
    df_rgb_conversion_table = pl.read_csv(
        rgb_conversion_table_file,
        separator="\t",
        comment_prefix="#",
        infer_schema_length=10000,
    )

    tec_df, converted_rgd_count = convert_rgd(tec_df, df_rgb_conversion_table)

    del df_rgb_conversion_table

    ## Process ens
    df_ens_conversion_table = pl.read_csv(
        ens_conversion_table_file,
        separator="\t",
        comment_prefix="#",
        infer_schema_length=10000,
    )

    tec_df, converted_ens_count = convert_ens(tec_df, df_ens_conversion_table)

    del df_ens_conversion_table


    ## Remove genes with multiple names
    is_unique_list: list[bool] = [
        is_unique(gene)
        for gene in tec_df["gene_symbol"]
    ]

    tec_df = tec_df.with_columns(
        pl.Series("is_unique", is_unique_list)
    )

    ### Collapse duplicates
    tec_df, dropped_duplicates_count = collapse_duplicate_genes_by_mean(tec_df, 0.1)
    result.dropped_duplicates_count = dropped_duplicates_count

    tec_df, unresolved_identifier_count = remove_unsolved_genes(tec_df)
    result.unresolved_identifiers = unresolved_identifier_count

    result.cleaned_genes_count = tec_df.height

    ##########################################################################################
    #                                       SAVE
    ##########################################################################################

    cleaned_output_file: Path = output / dataset.file.name
    tec_df.write_csv(cleaned_output_file)

    save_preprocessing_result(
        result=result,
        output_path=output,
        filename=dataset.file.stem,
    )
print("Done!")