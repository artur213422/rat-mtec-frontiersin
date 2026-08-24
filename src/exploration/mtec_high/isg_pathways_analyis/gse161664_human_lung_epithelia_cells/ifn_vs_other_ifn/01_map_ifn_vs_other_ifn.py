from pathlib import Path
from typing import Any, TypeAlias

import polars as pl

from src.common.paths import find_file_in_project, resolve_dir_from_root
from src.exploration.isg_pathways_analyis.gse161664_human_lung_epithelia_cells.internial.pick_isg import is_candidate, \
    CandidateResult

ResultValue: TypeAlias = str | float
ResultRow: TypeAlias = dict[str, ResultValue]

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

output_dir: Path = resolve_dir_from_root("data/analysis/isg_pathways_analysis/gse161664_human_lung_epithelia_cells/ifn_vs_other_ifn/isg_pathways")

output_dir.mkdir(parents=True, exist_ok=True)

output_file: Path = output_dir / "mapped_ifn_vs_other_ifn.csv"

input_file: Path = find_file_in_project("gse161664_counts.csv")
gse161664_df: pl.DataFrame = pl.read_csv(input_file)

IFN_NAMES: tuple[str, ...] = (
    "ifna",
    "ifnb",
    "ifng",
    "ifnl",
)


##########################################################################################
#                                       PROCESS DATA
##########################################################################################

result_rows: list[ResultRow] = []

for row in gse161664_df.iter_rows(named=True):
    typed_row: dict[str, Any] = row
    gene: str = str(typed_row["gene"])

    candidate: CandidateResult | None = is_candidate(row=typed_row, ifns=list(IFN_NAMES))

    if candidate is None:
        continue

    result_rows.append(
        {
            "gene": candidate.gene_symbol,
            "ifn": candidate.ifn,
            "log2_fold_change": candidate.log2_fold_change,
            "p_value": candidate.p_value,
        }
    )


result_df: pl.DataFrame = pl.DataFrame(
    result_rows,
    schema={
        "gene": pl.String,
        "ifn": pl.String,
        "log2_fold_change": pl.Float64,
        "p_value": pl.Float64,
    },
    orient="row",
).sort(
    by=["ifn", "gene"],
)


##########################################################################################
#                                           SAVE
##########################################################################################

result_df.write_csv(output_file)

print(f"Saved {result_df.height} IFN-specific genes to:")
