from pathlib import Path
from statistics import mean

import polars as pl

from src.common.dataframe_ops import get_genotype_expression_values
from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.paths import try_find_file_in_project, resolve_dir_from_root
from src.common.types.enums import Genotypes

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

out_file: Path = resolve_dir_from_root("data/mtec_high/analysis/tsa_expression_pattern") / "tsa_mtec_expression.csv"


ranked_tsa_file: Path = try_find_file_in_project("ranked_tsa.csv")
tsa_df: pl.DataFrame = pl.read_csv(ranked_tsa_file)

datasets: list[DatasetMetadata] = get_mtec_high()

dataset_frames: dict[str, pl.DataFrame] = {
    dataset.dataset_id: dataset.get_standardized_orthologs_df()
    for dataset in datasets
}
##########################################################################################
#                                       SETTINGS
##########################################################################################

top_tsa_count: int = 60
minimum_dataset_count: int = 2
required_dataset: str = "rat_mtec_high"

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

result_rows: list[dict[str, object]] = []

tissue_groups: list[str] = (
    tsa_df
    .get_column("tissue_group")
    .unique()
    .sort()
    .to_list()
)

for tissue_group in tissue_groups:
    print("Analyzing", tissue_group)

    group_tsa: pl.DataFrame = tsa_df.filter(pl.col("tissue_group") == tissue_group).sort("rank")

    selected_gene_count: int = 0

    for tsa_row in group_tsa.iter_rows(named=True):
        if selected_gene_count >= top_tsa_count:
            break

        gene_symbol: str = str(tsa_row["gene_symbol"])
        ortholog_id: str = str(tsa_row["ortholog_id"])

        matches_by_dataset: dict[str, pl.DataFrame] = {
            dataset_name: dataset_df.filter(
                pl.col("ortholog_id") == ortholog_id
            )
            for dataset_name, dataset_df in dataset_frames.items()
        }

        rat_matches: pl.DataFrame = matches_by_dataset[required_dataset]

        if rat_matches.height == 0:
            continue

        found_dataset_count: int = sum(
            matches.height == 1
            for matches in matches_by_dataset.values()
        )

        if found_dataset_count < minimum_dataset_count:
            continue

        selected_gene_count += 1
        rank: int = selected_gene_count

        for dataset_name, matches in matches_by_dataset.items():
            if matches.height > 1:
                raise ValueError(f"Expected at most one row for ortholog_id={ortholog_id}, found {matches.height} in dataset {dataset_name}")

            if matches.height == 0:
                result_rows.append(
                    {
                        "tissue_group": tissue_group,
                        "rank": rank,
                        "gene_symbol": gene_symbol,
                        "ortholog_id": ortholog_id,
                        "dataset": dataset_name,
                        "found_in_dataset": False,
                        "ctrl_mean": None,
                        "ko_mean": None,
                        "log2_fold_change": None,
                        "p_adj": None,
                    }
                )
                continue

            row: dict[str, object] = matches.row(0, named=True)

            ctrl_values: list[float] = get_genotype_expression_values(
                row,
                Genotypes.CTRL,
            )

            ko_values: list[float] = get_genotype_expression_values(
                row,
                Genotypes.KO,
            )

            result_rows.append(
                {
                    "tissue_group": tissue_group,
                    "rank": rank,
                    "gene_symbol": gene_symbol,
                    "ortholog_id": ortholog_id,
                    "dataset": dataset_name,
                    "found_in_dataset": True,
                    "ctrl_mean": mean(ctrl_values),
                    "ko_mean": mean(ko_values),
                    "log2_fold_change": float(row["log2_fold_change"]),
                    "p_adj": float(row["p_adj"]),
                }
            )

    if selected_gene_count < top_tsa_count:
        print(f"{tissue_group}: only {selected_gene_count} TSAs were found in at least {minimum_dataset_count} datasets")

result_df: pl.DataFrame = pl.DataFrame(result_rows)

##########################################################################################
#                                           SAVE
##########################################################################################

result_df.write_csv(out_file)