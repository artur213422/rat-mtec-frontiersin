from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import polars as pl

from config.paths import RSCRIPT_EXE
from src.common.paths import find_file_in_project, resolve_dir_from_root
from src.data_preprocessing.mtec_high.internal.R.run_R import run_rscript


LimmaDataType = Literal["rma", "rpkm"]


LIMMA_RESULT_COLUMNS: tuple[str, ...] = (
    "log2_fold_change",
    "t_statistic",
    "ranking_statistic",
    "p_val",
    "p_adj",
    "average_log2_expression",
    "test",
)


@dataclass(slots=True)
class LimmaInput:
    frame: pl.DataFrame
    output_file: Path
    data_type: LimmaDataType

    merge: bool = True

    rscript_exe: Path = RSCRIPT_EXE
    r_script: Path = find_file_in_project("limma_script.R")
    temporary_folder: Path = resolve_dir_from_root("data/temp/limma")

    pseudocount: float = 0.1

    def __post_init__(self) -> None:
        if not self.rscript_exe.is_file():
            raise FileNotFoundError(f"Rscript executable was not found: {self.rscript_exe}")

        if not self.r_script.is_file():
            raise FileNotFoundError("Limma R script was not found: {self.r_script}")

        if self.data_type == "rpkm" and self.pseudocount <= 0:
            raise ValueError("pseudocount must be positive for RPKM data")

        if "gene_symbol" not in self.frame.columns:
            raise ValueError("Input frame is missing 'gene_symbol'")


def run_limma(limma_input: LimmaInput) -> pl.DataFrame:
    control_columns: list[str] = sorted(
        column
        for column in limma_input.frame.columns
        if column.startswith("ctrl")
        and column.removeprefix("ctrl").isdigit()
    )

    ko_columns: list[str] = sorted(
        column
        for column in limma_input.frame.columns
        if column.startswith("ko")
        and column.removeprefix("ko").isdigit()
    )

    if len(control_columns) < 2:
        raise ValueError(f"At least two control columns are required; found {control_columns}")

    if len(ko_columns) < 2:
        raise ValueError(f"At least two KO columns are required; found {ko_columns}")


    temporary_input_file: Path = limma_input.temporary_folder / f"{limma_input.output_file.stem}_limma_input.csv"


    temporary_result_file: Path = limma_input.temporary_folder / f"{limma_input.output_file.stem}_limma_result.csv"


    limma_input.frame.select(
        "gene_symbol",
        *control_columns,
        *ko_columns,
    ).write_csv(
        temporary_input_file
    )

    config: dict[str, object] = {
        "input_file": str(temporary_input_file),
        "output_file": str(temporary_result_file),
        "data_type": limma_input.data_type,
        "pseudocount": limma_input.pseudocount,
        "control_columns": control_columns,
        "ko_columns": ko_columns,
    }

    run_rscript(
        rscript_exe=limma_input.rscript_exe,
        r_script=limma_input.r_script,
        config=config,
    )

    if not temporary_result_file.is_file():
        raise RuntimeError(f"Limma did not create the expected result file: {temporary_result_file}")

    # These are the columns actually written by the R script.
    limma_df: pl.DataFrame = pl.read_csv(
        temporary_result_file,
        schema_overrides={
            "gene_symbol": pl.String,
            "log2_fold_change": pl.Float64,
            "t_statistic": pl.Float64,
            "p_val": pl.Float64,
            "p_adj": pl.Float64,
            "average_log2_expression": pl.Float64,
        },
    )

    required_r_columns: set[str] = {
        "gene_symbol",
        "log2_fold_change",
        "t_statistic",
        "p_val",
        "p_adj",
        "average_log2_expression",
    }

    missing_r_columns: set[str] = (
        required_r_columns
        - set(limma_df.columns)
    )

    if missing_r_columns:
        raise ValueError(f"Limma result is missing columns: {sorted(missing_r_columns)}")

    # Create common cross-method columns in Python.
    limma_df = limma_df.with_columns(
        pl.col("t_statistic")
        .alias("ranking_statistic"),

        pl.lit("limma_moderated_t")
        .alias("test"),
    )

    if not limma_input.merge:
        result_df: pl.DataFrame = limma_df

    else:
        duplicate_genes_df: pl.DataFrame = (
            limma_input.frame
            .group_by("gene_symbol")
            .len()
            .filter(
                pl.col("len") > 1
            )
        )

        if duplicate_genes_df.height > 0:
            raise ValueError("Cannot merge limma results because the input contains duplicate gene symbols.")

        old_result_columns: list[str] = [
            column
            for column in LIMMA_RESULT_COLUMNS
            if column in limma_input.frame.columns
        ]

        unchanged_df: pl.DataFrame = (
            limma_input.frame.drop(
                old_result_columns
            )
            if old_result_columns
            else limma_input.frame
        )

        result_df = unchanged_df.join(
            limma_df.select(
                "gene_symbol",
                *LIMMA_RESULT_COLUMNS,
            ),
            on="gene_symbol",
            how="left",
            validate="1:1",
        )

        missing_statistics: int = (
            result_df
            .get_column("ranking_statistic")
            .null_count()
        )

        if missing_statistics > 0:
            raise ValueError(f"Limma results were missing for {missing_statistics:,} input genes.")

    result_df.write_csv(limma_input.output_file)

    return result_df