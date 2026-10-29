from pathlib import Path
import polars as pl
import pandas as pd
from src.common.paths import find_file_in_project, PROJECT_ROOT

excel_file: Path = find_file_in_project("NIHMS1834913-supplement-2.xlsx")

out_csv: Path = PROJECT_ROOT / "resources/michelson_2022_et_al_mimetic_cells/source/mimetic_cells.csv"

ignore_sheets: set[str] = {"Coversheet"}


def read_sheet(workbook: pd.ExcelFile, sheet_name: str) -> pl.DataFrame:
    df: pd.DataFrame = pd.read_excel(workbook, sheet_name=sheet_name)

    gene_column: str = next(
        column
        for column in df.columns
        if str(column).startswith("rownames(")
    )

    df = df.rename(columns={gene_column: "gene_symbol"})

    return (
        pl.from_pandas(df)
        .with_columns(
            pl.lit(sheet_name).alias("population")
        )
        .select(
            "population",
            "gene_symbol",
            pl.exclude("population", "gene_symbol"),
        )
    )


with pd.ExcelFile(excel_file) as workbook:
    sheet_names: list[str] = [
        sheet_name
        for sheet_name in workbook.sheet_names
        if sheet_name not in ignore_sheets
    ]

    frames: list[pl.DataFrame] = [
        read_sheet(workbook, sheet_name)
        for sheet_name in sheet_names
    ]



mimetic_cells: pl.DataFrame = pl.concat(frames, how="diagonal_relaxed")
mimetic_cells = mimetic_cells.rename({
    "avg_log2FC": "log2_fold_change",
    "p_val_adj": "p_adj",
    "pct.1": "pct_1",
    "pct.2": "pct_2"
})

mimetic_cells.write_csv(out_csv)
