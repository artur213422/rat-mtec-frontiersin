import re
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.table import TableStyleInfo, Table


def make_excel_table_name(name: str) -> str:
    table_name = re.sub(r"[^A-Za-z0-9_]", "_", name)

    if not table_name:
        table_name = "Table1"

    if table_name[0].isdigit():
        table_name = f"T_{table_name}"

    return table_name[:255]

def create_excel_table(excel_path: Path, sheet_name: str = "Sheet1",  table_name: str | None = None) -> Path:
    if not excel_path.exists():
        raise FileNotFoundError(excel_path)

    wb = load_workbook(excel_path)

    if sheet_name not in wb.sheetnames:
        raise ValueError(
            f"Sheet '{sheet_name}' not found. Available sheets: {wb.sheetnames}"
        )

    ws = wb[sheet_name]

    if ws.max_row < 2:
        wb.save(excel_path)
        print(f"Saved Excel file without table because sheet has no data rows: {excel_path}")
        return excel_path

    if ws.max_column < 1:
        wb.save(excel_path)
        print(f"Saved Excel file without table because sheet has no columns: {excel_path}")
        return excel_path

    # Excel table headers must be non-empty strings.
    for column_idx in range(1, ws.max_column + 1):
        cell = ws.cell(row=1, column=column_idx)

        if cell.value is None or str(cell.value).strip() == "":
            cell.value = f"Column_{column_idx}"
        else:
            cell.value = str(cell.value)

    last_col = get_column_letter(ws.max_column)
    table_ref = f"A1:{last_col}{ws.max_row}"

    if table_name is None:
        table_name = make_excel_table_name(excel_path.stem)

    table = Table(
        displayName=table_name,
        ref=table_ref,
    )

    style = TableStyleInfo(
        name="TableStyleMedium9",
        showFirstColumn=False,
        showLastColumn=False,
        showRowStripes=True,
        showColumnStripes=False,
    )

    table.tableStyleInfo = style
    ws.add_table(table)

    ws.freeze_panes = "A2"

    for column_idx in range(1, ws.max_column + 1):
        column_letter = get_column_letter(column_idx)

        values = [
            ws.cell(row=row_idx, column=column_idx).value
            for row_idx in range(1, min(ws.max_row, 101) + 1)
        ]

        max_length = max(len(str(value)) for value in values if value is not None)

        ws.column_dimensions[column_letter].width = min(
            max(max_length + 2, 10),
            40,
        )

    wb.save(excel_path)

    print(f"Saved Excel file with table: {excel_path}")

    return excel_path


def make_excel_table_name(name: str) -> str:
    table_name = re.sub(r"[^A-Za-z0-9_]", "_", name)

    if not table_name:
        table_name = "Table1"

    if table_name[0].isdigit():
        table_name = f"T_{table_name}"

    return table_name[:255]