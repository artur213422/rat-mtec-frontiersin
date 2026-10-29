from pathlib import Path
import re

import polars as pl
from src.common.feratures.limma_voom.run_limma_voom import LimmaVoomInput, run_limma_voom
from src.common.paths import resolve_dir_from_root, find_file_in_project

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

source_count_file: Path = find_file_in_project("gse161664_counts.csv")
source_count_df: pl.DataFrame = pl.read_csv(source_count_file)

prepared_for_limma_folder: Path = resolve_dir_from_root("resources/interferon/gse161664_lung/intermediate/ifn_vs_control/prepared_for_limma")

limma_result_folder: Path = resolve_dir_from_root("resources/interferon/gse161664_lung/intermediate/ifn_vs_control/limma_result")

SAMPLE_PREFIXES = {
    "ctrl": "ctrl",
    "ifna": "ifna",
    "ifnb": "ifnb",
    "ifng": "ifng",
    "ifnl": "ifnl3",
}

IFN_NAMES = ("ifna", "ifnb", "ifng", "ifnl")


##########################################################################################
#                                       PROCESS DATA
##########################################################################################
def sample_columns(condition: str) -> list[str]:
    prefix = SAMPLE_PREFIXES[condition]

    columns = [
        column
        for column in source_count_df.columns
        if re.fullmatch(rf"{re.escape(prefix)}\d+", column)
    ]

    if len(columns) < 2:
        raise ValueError(f"Expected at least two samples for {condition}; found {columns}")

    return sorted(
        columns,
        key=lambda column: int(column[len(prefix):]),
    )



control_columns = sample_columns("ctrl")

for ifn_name in IFN_NAMES:
    target_columns = sample_columns(ifn_name)
    selected_columns = control_columns + target_columns

    # Convert wide counts to one row per gene/sample combination.
    counts_df = source_count_df.select(
        "gene_symbol",
        *selected_columns,
    ).unpivot(
        index="gene_symbol",
        on=selected_columns,
        variable_name="sample_id",
        value_name="count",
    )

    # Row order matches selected_columns.
    metadata_df = pl.DataFrame({
        "sample_id": selected_columns,
        "experiment": [1] * len(selected_columns),
        "condition": (
            ["ctrl"] * len(control_columns)
            + [ifn_name] * len(target_columns)
        ),
    })

    counts_file = prepared_for_limma_folder / f"{ifn_name}_vs_ctrl_counts.csv"

    metadata_file = prepared_for_limma_folder / f"{ifn_name}_vs_ctrl_metadata.csv"

    output_file = limma_result_folder / f"{ifn_name}_vs_ctrl.csv"


    counts_df.write_csv(counts_file)
    metadata_df.write_csv(metadata_file)

    run_limma_voom(
        LimmaVoomInput(
            counts_file=counts_file,
            metadata_file=metadata_file,
            output_file=output_file,
            reference_condition="ctrl",
            target_condition=ifn_name,
        )
    )

    print(f"Saved {ifn_name} versus control: {output_file}")