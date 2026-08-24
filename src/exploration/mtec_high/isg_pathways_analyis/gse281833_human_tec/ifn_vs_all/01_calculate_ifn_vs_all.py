from pathlib import Path
from typing import Any, TypeAlias

import polars as pl

from src.common.dataframe_ops import get_custom_genotype_expression_values, get_genotype_expression_values
from src.common.feratures.math.fold_change import calculate_log2_fold_change
from src.common.feratures.math.is_significant import is_significant
from src.common.feratures.math.logarithm import log2_transform
from src.common.feratures.math.p_val import calculate_raw_p_val
from src.common.paths import find_file_in_project, resolve_dir_from_root
from src.common.types.enums import Genotypes


ResultValue: TypeAlias = str | float
ResultRow: TypeAlias = dict[str, ResultValue]

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

output_dir: Path = resolve_dir_from_root("data/analysis/isg_pathways_analysis/gse281833_human_tec/ifn_vs_all_genes/isg_pathways")

output_file: Path = output_dir / "gse281833_mapped_ifn_vs_all.csv"

input_file: Path = find_file_in_project("gse281833_counts.csv")
gse161664_df: pl.DataFrame = pl.read_csv(input_file)

IFN_NAMES: tuple[str, ...] = ("ifna", "ifng")


##########################################################################################
#                                       PROCESS DATA
##########################################################################################

result_rows: list[ResultRow] = []

for row in gse161664_df.iter_rows(named=True):
    typed_row: dict[str, Any] = row
    gene: str = str(typed_row["gene"])

    raw_ctrl_values: list[float] = get_genotype_expression_values(row=typed_row, gt=Genotypes.CTRL)
    log2_ctrl_values: list[float] = log2_transform(raw_ctrl_values)

    for ifn_name in IFN_NAMES:
        raw_ifn_values: list[float] = (get_custom_genotype_expression_values(row=typed_row, pattern=ifn_name))

        log2_ifn_values: list[float] = log2_transform(raw_ifn_values)

        log2_fc: float = calculate_log2_fold_change(ctrl_values=log2_ctrl_values, ko_values=log2_ifn_values)
        p_value: float = calculate_raw_p_val(ctrl_values=log2_ctrl_values, ko_values=log2_ifn_values)

        if not is_significant(log2_fold_change=log2_fc, p_value=p_value):
            continue

        if log2_fc <= 0.0:
            continue

        result_rows.append(
            {
                "gene": gene,
                "ifn": ifn_name,
                "log2_fold_change": log2_fc,
                "p_value": p_value,
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
    by=["ifn", "log2_fold_change"],
    descending=[False, True],
)


##########################################################################################
#                                           SAVE
##########################################################################################

result_df.write_csv(output_file)

print(result_df)
print(f"Saved {result_df.height} gene–IFN associations to:")
print(output_file)