from pathlib import Path

import polars as pl

from src.common.excel import create_excel_table
from src.common.file_io import save_csv_to_excel
from src.common.paths import PROJECT_ROOT

input_dir =  PROJECT_ROOT / "data/mtec_high/derived_data/gsea_tsa_excluded/gsea/standard"
output =  PROJECT_ROOT / "data/mtec_high/derived_data/gsea_tsa_excluded/gsea/standard/selected"


input_files = (
    "go_biological_process_gsea.csv",
    "go_molecular_function_gsea.csv",
    "go_cellular_component_gsea.csv",
)

for input_name in input_files:
    input_file = input_dir / input_name
    gsea_df = pl.read_csv(input_file)

    retained_terms = (
        gsea_df
        .filter(pl.col("nominal_p_value") < 0.05)
        .select("term")
        .unique()
    )

    # Keep every dataset row for a pathway if any dataset has p < 0.05.
    filtered_df = (
        gsea_df
        .join(retained_terms, on="term", how="semi")
        .sort(["term", "dataset_id"])
    )

    output_file = output / input_name.replace(
        ".csv",
        "_one_dataset_nominal_p_lt_0_05.csv",
    )

    filtered_df.write_csv(output_file)

    output_excel: Path = save_csv_to_excel(output_file)

    create_excel_table(output_excel)

    print(f"{input_name}: retained {retained_terms.height} pathways; saved {filtered_df.height} rows to {output_file}")