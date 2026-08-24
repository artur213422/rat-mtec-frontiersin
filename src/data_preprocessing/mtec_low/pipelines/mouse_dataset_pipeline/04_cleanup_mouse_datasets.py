from pathlib import Path

import polars as pl

from src.common.dataset_catalog.api.get_datasets import get_mouse_tec_low
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.paths import resolve_dir_from_root, find_file_in_project
from src.data_preprocessing.mtec_high.internal.drop_low_expressing_genes import drop_low_expressing_genes
from src.data_preprocessing.mtec_high.internal.gene_symbol_cleaning_tools import unify_gene_symbol_splitters, drop_invalid_values, \
    cleanup_gene_symbol, cleanup_mouse_rik_genes, collapse_duplicate_genes_by_mean, remove_unsolved_genes
from src.data_preprocessing.mtec_high.internal.models.preprocessing_result.model import PreprocessingResult
from src.data_preprocessing.mtec_high.internal.models.preprocessing_result.ops import save_preprocessing_result
from src.data_preprocessing.mtec_high.internal.rat_ens_and_rgb_conversion import convert_ens

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

mouse_datasets: list[DatasetMetadata] = get_mouse_tec_low()


source_folder: Path = resolve_dir_from_root("data/mtec_low/raw/intermediate/annotated_expression_matrix/mouse")
output_folder: Path = resolve_dir_from_root("data/mtec_low/raw/intermediate/cleaned_names_expression_matrix/mouse")

rik_annotation_table_file: Path = find_file_in_project("MRK_List1.rpt")
ens_conversion_table_file: Path = find_file_in_project("biomart_ens_mouse_symbol_annotations.csv")

##########################################################################################
#                                       PROCESS FILE
##########################################################################################
for dataset in mouse_datasets:
    ##  Preparation
    result: PreprocessingResult = PreprocessingResult(dataset.dataset_name)


    source_file: Path = source_folder / f"{dataset.dataset_id}.csv"
    print(f"Running {source_file.stem}")

    if not source_file.exists():
        raise FileNotFoundError(f"{source_file} does not exist")

    tec_df: pl.DataFrame = pl.read_csv(source_file)

    result.total_gene_count = tec_df.height

    ### Unify
    tec_df = unify_gene_symbol_splitters(tec_df)

    #tec_df = standardize_sample_group_col_names(tec_df, dataset.ctrl_cols, dataset.ko_cols)

    ### Drop nulls
    tec_df, dropped_invalid_count = drop_invalid_values(tec_df)

    result.dropped_invalid_count = dropped_invalid_count

    ## Drop zero expression genes

    tec_df = drop_low_expressing_genes(tec_df, 0.0)

    ## Process ens
    df_ens_conversion_table = pl.read_csv(
        ens_conversion_table_file,
        separator=",",
        infer_schema_length=10000,
    )

    tec_df, converted_ens_count = convert_ens(tec_df, df_ens_conversion_table)

    del df_ens_conversion_table

    ## Clean symbols
    cleaned_symbols: list[str] = []

    for gene in tec_df["gene_symbol"]:
        cleaned_gene: str = cleanup_gene_symbol(gene)
        cleaned_symbols.append(cleaned_gene)

    tec_df = tec_df.with_columns(
        pl.Series("gene_symbol", cleaned_symbols),
    )


    tec_df, resolved_rik_count = cleanup_mouse_rik_genes(tec_df, rik_annotation_table_file)

    ### Collapse duplicates
    tec_df, dropped_duplicates_count = collapse_duplicate_genes_by_mean(tec_df, 0.1)
    result.dropped_duplicates_count = dropped_duplicates_count

    tec_df, unresolved_identifier_count = remove_unsolved_genes(tec_df)
    result.unresolved_identifiers = unresolved_identifier_count

    result.cleaned_genes_count = tec_df.height


    ##########################################################################################
    #                                       SAVE
    ##########################################################################################

    cleaned_output_file: Path = output_folder / f"{dataset.dataset_id}.csv"
    tec_df.write_csv(cleaned_output_file)

    save_preprocessing_result(
        result=result,
        output_path=output_folder,
        filename=dataset.dataset_id,
    )


