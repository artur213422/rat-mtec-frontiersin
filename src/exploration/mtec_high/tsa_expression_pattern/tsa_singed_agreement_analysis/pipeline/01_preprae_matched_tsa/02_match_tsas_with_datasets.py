from pathlib import Path

import polars as pl

from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.paths import try_find_file_in_project, resolve_dir_from_root

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

defined_tsa_file: Path = try_find_file_in_project("defined_tsa.csv")
defined_tsa_df: pl.DataFrame = pl.read_csv(defined_tsa_file)

defined_tsa_df = defined_tsa_df.with_columns(
    pl.col("ortholog_id").cast(pl.String)
)

matched_tsa_file: Path = resolve_dir_from_root("data/mtec_high/analysis/tsa_expression_pattern/tsa_singed_agreement_analysis/tsa_list/intermediate") / "matched_tsa.csv"

datasets: list[DatasetMetadata] = get_mtec_high()

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

result_rows: list[dict[str, object]] = []

for dataset in datasets:
    dataset_name: str = dataset.dataset_id
    mtec_df: pl.DataFrame = dataset.get_standardized_df()
    print("Analysing dataset: ", dataset_name)
    for tsa_row in defined_tsa_df.iter_rows(named=True):
        gene_symbol: str = str(tsa_row["gene_symbol"])
        ortholog_id: str | None = tsa_row["ortholog_id"]
        tissue_group: str = str(tsa_row["tissue_group"])
        rank: int = int(tsa_row["rank"])

        matches: pl.DataFrame = pl.DataFrame()
        matched_by: str | None = None


        # First try ortholog ID.
        if ortholog_id is not None:
            matches = mtec_df.filter(
                (pl.col("ortholog_id") == ortholog_id) & (pl.col("has_clean_ortholog") == True)
            )

            if matches.height > 0:
                matched_by = "ortholog_id"

        # Fall back to gene symbol.
        if matches.height == 0:
            matches = mtec_df.filter(
                pl.col("gene_symbol") == gene_symbol
            )

            if matches.height > 0:
                matched_by = "gene_symbol"

        if matches.height > 1:
            raise ValueError(f"Multiple matches for TSA '{gene_symbol}' in dataset '{dataset_name}'" )

        if matches.height == 0:
            result_rows.append(
                {
                    "dataset": dataset_name,
                    "tissue_group": tissue_group,
                    "rank": rank,
                    "gene_symbol": gene_symbol,
                    "ortholog_id": ortholog_id,
                    "found_in_dataset": False,
                    "matched_by": None,
                }
            )
            continue

        dataset_row: dict[str, object] = matches.row(
            0,
            named=True,
        )

        # TSA identity should stay authoritative.
        dataset_row.pop("gene_symbol", None)
        dataset_row.pop("ortholog_id", None)

        result_rows.append(
            {
                "dataset": dataset_name,
                "tissue_group": tissue_group,
                "rank": rank,
                "gene_symbol": gene_symbol,
                "ortholog_id": ortholog_id,
                "found_in_dataset": True,
                "matched_by": matched_by,
                **dataset_row,
            }
        )

matched_tsa_df: pl.DataFrame = pl.DataFrame(result_rows)

##########################################################################################
#                                           SAVE
##########################################################################################

matched_tsa_df.write_csv(matched_tsa_file)