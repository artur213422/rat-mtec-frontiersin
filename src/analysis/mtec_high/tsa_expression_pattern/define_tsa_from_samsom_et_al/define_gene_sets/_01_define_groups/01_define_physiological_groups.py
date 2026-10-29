import json
from pathlib import Path
import sys

# Support direct execution from any working directory.
PROJECT_ROOT = Path(__file__).resolve().parents[7]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.analysis.mtec_high.tsa_expression_pattern.define_tsa_from_samsom_et_al.define_gene_sets.config.paths import (
    ATLAS_FILE,
    TISSUE_GROUP_FILE,
    GROUPED_ATLAS_FILE,
)

import polars as pl


##########################################################################################
#                                       PREPARE DATA
##########################################################################################

atlas_file: Path = ATLAS_FILE
atlas_df: pl.DataFrame = pl.read_csv(atlas_file)

tissue_group_file: Path = TISSUE_GROUP_FILE

out_file: Path = GROUPED_ATLAS_FILE
out_file.parent.mkdir(parents=True, exist_ok=True)

##########################################################################################
# PROCESS DATA
##########################################################################################

metadata_columns: list[str] = [
    "ID_REF",
    "ortholog_id",
    "gene_symbol",
    "has_clean_ortholog",
]

with tissue_group_file.open("r", encoding="utf-8",) as file:
    tissue_groups: dict[str, list[str]] = json.load(file)

if not tissue_groups:
    raise ValueError("tissue_groups.json contains no tissue groups")

empty_groups: list[str] = [
    group_name
    for group_name, member_tissues in tissue_groups.items()
    if not member_tissues
]

if empty_groups:
    raise ValueError(f"Tissue groups without members: {empty_groups}")

all_member_tissues: list[str] = [
    tissue
    for member_tissues in tissue_groups.values()
    for tissue in member_tissues
]

duplicated_tissues: list[str] = sorted(
    {
        tissue
        for tissue in all_member_tissues
        if all_member_tissues.count(tissue) > 1
    }
)

if duplicated_tissues:
    raise ValueError(f"Tissues assigned to multiple groups: {duplicated_tissues}")

missing_tissues: list[str] = sorted(set(all_member_tissues)  - set(atlas_df.columns))

if missing_tissues:
    raise ValueError(f"Tissues missing from the atlas: {missing_tissues}")

group_expressions: list[pl.Expr] = []

for group_name, member_tissues in tissue_groups.items():
    # Convert log2 RMA to linear expression.
    linear_member_expressions: list[pl.Expr] = [
        pl.lit(2.0).pow(
            pl.col(tissue)
        )
        for tissue in member_tissues
    ]

    # Calculate arithmetic mean on the linear scale,
    # then return it to log2.
    group_expression: pl.Expr = (
        pl.mean_horizontal(
            *linear_member_expressions
        )
        .log(2.0)
        .alias(group_name)
    )

    group_expressions.append(group_expression)

atlas_df = atlas_df.select(*metadata_columns, *group_expressions,)

print(f"Final tissue groups: {len(tissue_groups)}")

##########################################################################################
# SAVE
##########################################################################################

atlas_df.write_csv(out_file)

print(f"Saved: {out_file}")