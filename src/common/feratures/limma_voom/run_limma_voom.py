import polars as pl

from dataclasses import dataclass
from pathlib import Path
from typing import Literal
from config.paths import RSCRIPT_EXE
from src.common.paths import find_file_in_project
from src.data_preprocessing.mtec_high.internal.R.run_R import run_rscript

LIMMA_VOOM_REQUIRED_COUNTS_COLUMNS = frozenset({
    "gene_symbol",
    "sample_id",
    "count",
})

LIMMA_VOOM_REQUIRED_METADATA_COLUMNS = frozenset({
    "sample_id",
    "experiment",
    "condition",
})


@dataclass(frozen=True, slots=True)
class LimmaVoomInput:
    counts_file: Path
    metadata_file: Path
    output_file: Path
    reference_condition: str # 'ctrl' etc
    target_condition: str # 'ko' etc
    replicate_type: Literal["biological"] = "biological"


    def __post_init__(self) -> None:
        if self.replicate_type != "biological":
            raise ValueError("LimmaVoomInput requires biological replicates. Technical replicates are not supported.")

        counts_df = pl.read_csv(self.counts_file)
        metadata_df = pl.read_csv(self.metadata_file)

        missing_counts_cols = LIMMA_VOOM_REQUIRED_COUNTS_COLUMNS - set(counts_df.columns)
        if missing_counts_cols:
            raise ValueError(f"Counts file is missing required columns: {sorted(missing_counts_cols)}")

        missing_metadata_cols = LIMMA_VOOM_REQUIRED_METADATA_COLUMNS - set(metadata_df.columns)
        if missing_metadata_cols:
            raise ValueError(f"Metadata file is missing required columns: {sorted(missing_metadata_cols)}")


    def to_config(self) -> dict[str, object]:
        return {
            "counts_file": str(self.counts_file),
            "metadata_file": str(self.metadata_file),
            "output_file": str(self.output_file),
            "reference_condition": self.reference_condition,
            "target_condition": self.target_condition,
        }




def run_limma_voom(limma_input: LimmaVoomInput) -> None:
    r_script: Path = find_file_in_project("limma_voom.R")

    run_rscript(
        rscript_exe=RSCRIPT_EXE,
        r_script=r_script,
        config=limma_input.to_config(),
    )
