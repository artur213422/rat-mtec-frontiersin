from pathlib import Path
import polars as pl
from src.common.paths import resolve_file_from_root, resolve_dir_from_root

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


source_file: Path = resolve_file_from_root("resources/interferon/gse28183_tec/intermediate/ifn_vs_all/gse281833_tec_cleaned.csv")

output_file: Path = resolve_dir_from_root("resources/interferon/gse28183_tec/standardized/candidates") / "gse281833_ifn_vs_control_list.csv"


##########################################################################################
#                                       PROCESS
##########################################################################################

result_df = (
    pl.read_csv(source_file, schema_overrides={"ortholog_id": pl.String})
    .drop_nulls(subset=["ortholog_id"])
    .with_columns(pl.col("ortholog_id").str.split("///"))
    .explode("ortholog_id")
)

output_file.parent.mkdir(parents=True, exist_ok=True)
result_df.write_csv(output_file)

print(f"Saved: {output_file}")