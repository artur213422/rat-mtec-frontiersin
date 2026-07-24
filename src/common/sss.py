from pathlib import Path
from typing import Any

import polars as pl

from src.common.dataframe_ops import (
    get_genotype_expression_values,
)
from src.common.excel import create_excel_table
from src.common.feratures.math.p_val import (
    calculate_log2_fold_change,
    calculate_raw_p_val,
)
from src.common.file_io import save_csv_to_excel
from src.common.models.dataset_metadata.model import (
    DatasetCatalog,
)
from src.common.models.dataset_metadata.ops import (
    build_dataset_catalog,
)
from src.common.paths import resolve_dir_from_root
from src.common.types.enums import Genotypes


LOG2_FC_COLUMN: str = "log2_fold_change"
RAW_P_VALUE_COLUMN: str = "raw_p_value"


def add_gene_statistics(
    tec_df: pl.DataFrame,
) -> pl.DataFrame:
    log2_fold_changes: list[float] = []
    raw_p_values: list[float] = []

    for row in tec_df.iter_rows(named=True):
        typed_row: dict[str, Any] = row

        ctrl_values: list[float] = (
            get_genotype_expression_values(
                row=typed_row,
                gt=Genotypes.CTRL,
            )
        )

        ko_values: list[float] = (
            get_genotype_expression_values(
                row=typed_row,
                gt=Genotypes.KO,
            )
        )

        log2_fold_change: float = (
            calculate_log2_fold_change(
                ctrl_values,
                ko_values,
            )
        )

        raw_p_value: float = calculate_raw_p_val(
            ctrl_values,
            ko_values,
        )

        log2_fold_changes.append(log2_fold_change)
        raw_p_values.append(raw_p_value)

    return tec_df.with_columns(
        pl.Series(
            name=LOG2_FC_COLUMN,
            values=log2_fold_changes,
            dtype=pl.Float64,
        ),
        pl.Series(
            name=RAW_P_VALUE_COLUMN,
            values=raw_p_values,
            dtype=pl.Float64,
        ),
    )


def save_dataset_statistics(
    dataset_id: str,
    tec_df: pl.DataFrame,
    output_folder: Path,
) -> tuple[Path, Path]:
    result_df: pl.DataFrame = add_gene_statistics(
        tec_df=tec_df,
    )

    csv_file: Path = (
        output_folder
        / f"{dataset_id}_gene_statistics.csv"
    )

    result_df.write_csv(csv_file)

    excel_file: Path = save_csv_to_excel(csv_file)
    create_excel_table(excel_file)

    return csv_file, excel_file


def main() -> None:
    output_folder: Path = resolve_dir_from_root(
        "data/analysis/general_statistics/data"
    )

    output_folder.mkdir(
        parents=True,
        exist_ok=True,
    )

    datasets: DatasetCatalog = build_dataset_catalog()

    for dataset in datasets.datasets:
        dataset_id: str = dataset.dataset_id

        print(f"Processing dataset: {dataset_id}")

        tec_df: pl.DataFrame = (
            dataset.get_standardized_df()
        )

        csv_file, excel_file = save_dataset_statistics(
            dataset_id=dataset_id,
            tec_df=tec_df,
            output_folder=output_folder,
        )

        print(f"Saved CSV:   {csv_file}")
        print(f"Saved Excel: {excel_file}")


if __name__ == "__main__":
    main()