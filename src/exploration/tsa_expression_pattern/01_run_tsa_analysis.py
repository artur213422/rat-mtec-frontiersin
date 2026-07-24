from pathlib import Path
from statistics import mean

import polars as pl

from src.common.dataframe_ops import get_genotype_expression_values
from src.common.feratures.math.p_val import calculate_log2_fold_change, calculate_raw_p_val
from src.common.models.dataset_metadata.model import DatasetCatalog
from src.common.models.dataset_metadata.ops import build_dataset_catalog
from src.common.paths import try_find_file_in_project, resolve_dir_from_root
from src.common.types.enums import Genotypes
from src.exploration.tsa_expression_pattern.internal.tsa_list import DefinedTsa
from src.resource_processing.define_tsa_pipeline.internal.tissue_extraction import \
    extract_list_from_json

##########################################################################################
#                                       PREPARE DATA
##########################################################################################




picked_tissues_file: Path = try_find_file_in_project("picked_tissues.json")
tissues: list[str] = extract_list_from_json(picked_tissues_file)
out_file: Path = resolve_dir_from_root("data/analysis") / "tsa_expression_pattern/analysed_tsa.csv"
datasets: DatasetCatalog = build_dataset_catalog()

##########################################################################################
#                                       PROCESS DATA
##########################################################################################
result_rows: list[dict[str, object]] = []

for tissue in tissues:
    file: Path = try_find_file_in_project(f"{tissue}_tsa_candidates.csv")

    tsa: DefinedTsa = DefinedTsa.from_file(tissue, file)

    if len(tsa.ranked_tsa) == 0:
        print(f"Skipping {tissue}")
        continue

    for rank, gene_ref in enumerate(tsa.ranked_tsa, start=1):
        for dataset in datasets.datasets:
            dataset_df: pl.DataFrame = dataset.get_standardized_orthologs_df()
            dataset_name: str = dataset.dataset_id

            matches: pl.DataFrame = dataset_df.filter(
                pl.col("ortholog_id") == gene_ref.ortholog_id
            )

            if matches.height == 0:
                result_rows.append(
                    {
                        "tissue": tissue,
                        "rank": rank,
                        "gene_symbol": gene_ref.gene_symbol,
                        "ortholog_id": gene_ref.ortholog_id,
                        "dataset": dataset_name,
                        "found_in_dataset": False,
                        "ctrl_mean": None,
                        "ko_mean": None,
                        "log2_fold_change": None,
                        "p_value": None,
                    }
                )
                continue

            if matches.height > 1:
                raise ValueError(
                    f"Expected at most one row for ortholog_id={gene_ref.ortholog_id}, "
                    f"found {matches.height} in dataset {dataset_name}"
                )

            row: dict[str, object] = matches.row(0, named=True)

            ctrl_values: list[float] = get_genotype_expression_values(row, Genotypes.CTRL)
            ko_values: list[float] = get_genotype_expression_values(row, Genotypes.KO)

            fc: float = calculate_log2_fold_change(ctrl_values, ko_values)
            p_value: float = calculate_raw_p_val(ctrl_values, ko_values)

            result_rows.append(
                {
                    "tissue": tissue,
                    "rank": rank,
                    "gene_symbol": gene_ref.gene_symbol,
                    "ortholog_id": gene_ref.ortholog_id,
                    "dataset": dataset_name,
                    "found_in_dataset": True,
                    "ctrl_mean": mean(ctrl_values),
                    "ko_mean": mean(ko_values),
                    "log2_fold_change": fc,
                    "p_value": p_value,
                }
            )


result_df: pl.DataFrame = pl.DataFrame(result_rows)



##########################################################################################
#                                           SAVE
##########################################################################################

result_df.write_csv(out_file)

