import csv
import re
from pathlib import Path

import polars as pl

from src.common.paths import find_file_in_project

##########################################################################################
# PREPARE DATA
##########################################################################################
input_file: Path = find_file_in_project("GSE10246_affy_rma_log2_probe_by_sample.csv")

data: pl.DataFrame = pl.read_csv(input_file)

##########################################################################################
# VALIDATE
##########################################################################################

errors: list[str] = []

# Validate the original CSV header before relying on Polars.
with input_file.open("r", encoding="utf-8", newline="") as file:
    header: list[str] = next(csv.reader(file))

duplicated_columns: list[str] = sorted(
    {
        column
        for column in header
        if header.count(column) > 1
    }
)

if duplicated_columns:
    errors.append(f"Duplicate column names: {duplicated_columns}")

# ID_REF must exist.
if "ID_REF" not in data.columns:
    errors.append("Required column is missing: ID_REF")

# ID_REF should be the first column.
elif data.columns[0] != "ID_REF":
    errors.append(f"ID_REF must be the first column, but found: {data.columns[0]}")

# Every expression column must be a GSM accession.
sample_columns: list[str] = [
    column
    for column in data.columns
    if column != "ID_REF"
]

invalid_sample_columns: list[str] = [
    column
    for column in sample_columns
    if re.fullmatch(r"GSM\d+", column) is None
]

if invalid_sample_columns:
    errors.append(f"Invalid GSM sample column names: {invalid_sample_columns}")

# No nulls are allowed anywhere.
null_counts: dict[str, int] = data.null_count().row(
    0,
    named=True,
)

columns_with_nulls: dict[str, int] = {
    column: count
    for column, count in null_counts.items()
    if count > 0
}

if columns_with_nulls:
    errors.append(f"Null values found: {columns_with_nulls}")

# Probe IDs must not be null, empty or duplicated.
if "ID_REF" in data.columns:
    empty_probe_count: int = data.filter(
        pl.col("ID_REF").cast(pl.String).str.strip_chars() == ""
    ).height

    if empty_probe_count > 0:
        errors.append(f"Found {empty_probe_count} empty ID_REF values")

    duplicated_probe_ids: list[str] = (
        data
        .filter(pl.col("ID_REF").is_duplicated())
        .select("ID_REF")
        .unique()
        .sort("ID_REF")
        .get_column("ID_REF")
        .to_list()
    )

    if duplicated_probe_ids:
        preview = duplicated_probe_ids[:20]

        errors.append(f"Found {len(duplicated_probe_ids)} duplicated ID_REF values. First examples: {preview}")

# Expression columns must be numeric.
non_numeric_columns: list[str] = [
    column
    for column in sample_columns
    if not data.schema[column].is_numeric()
]

if non_numeric_columns:
    errors.append(f"Non-numeric expression columns: {non_numeric_columns}")

##########################################################################################
# REPORT
##########################################################################################

if errors:
    print("VALIDATION FAILED!")

    for error in errors:
        print(f"- {error}")

    raise ValueError(f"RMA matrix failed validation with {len(errors)} error(s)")

print("OK")
print(f"Probes: {data.height}")
print(f"GSM samples: {len(sample_columns)}")
print(f"Shape: {data.shape}")