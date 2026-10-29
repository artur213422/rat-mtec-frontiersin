from pathlib import Path
from typing import Any
import polars as pl

from src.common.paths import resolve_dir_from_root, find_file_in_project, collect_files_in_folder
from src.data_preprocessing.mtec_high.internal.R.run_R import run_rscript
from src.data_preprocessing.mtec_high.pipelines.mouse_dataset_pipeline.interials.helpers import gsm_sample_name_from_cel

##########################################################################################
#                                       PREPARE DATA
##########################################################################################
raw_cel_folder: Path = resolve_dir_from_root("resources/mouse_gene_atlas/raw/GSE10246_RAW")
cel_files: list[Path] = collect_files_in_folder(raw_cel_folder)

txt_output: Path = resolve_dir_from_root("resources/mouse_gene_atlas/raw") / "GSE10246_affy_rma_log2_probe_by_sample.txt"
csv_output: Path = txt_output.with_suffix(".csv")

r_exe = Path("C:/R/R/bin/Rscript.exe")

affy_script: Path = find_file_in_project("cel_to_txt_affy.R")
##########################################################################################
#                                      SETTINGS
##########################################################################################
platform_package: str = "mouse4302cdf"

##########################################################################################
#                                       PROCESS DATA
##########################################################################################
sample_names: list[str] = [
    gsm_sample_name_from_cel(cel_file)
    for cel_file in cel_files
]


if len(sample_names) != len(set(sample_names)):
    duplicated_names = sorted(
        {
            sample_name
            for sample_name in sample_names
            if sample_names.count(sample_name) > 1
        }
    )
    raise ValueError(f"Duplicate GSM accessions extracted from CEL filenames: {duplicated_names}")

config: dict[str, Any] = {
    "cel_files": [str(path) for path in cel_files],
    "sample_names": sample_names,
    "output_file": str(txt_output),
    "platform_package": platform_package
    ,"quote_output": True

}

run_rscript(r_exe, affy_script, config)
##########################################################################################
#                                       SAVE AS CSV
##########################################################################################


# The R script writes a tab-separated table.
matrix: pl.DataFrame = pl.read_csv(txt_output, separator="\t")

# Write a genuine comma-separated CSV.
matrix.write_csv(csv_output)

if not csv_output.is_file():
    raise RuntimeError(f"CSV output was not created: {csv_output}")

# Delete the temporary TXT only after successful conversion.
txt_output.unlink()

print(f"Saved: {csv_output}")
print(f"Deleted temporary file: {txt_output}")