from pathlib import Path

import polars as pl

from dataclasses import dataclass
from typing import Literal, Sequence
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
class MultiConditionalLimmaVoomInput:
    counts_file: Path
    metadata_file: Path
    output_file: Path
    contrasts_file: Path
    induction_log2_fc: float
    specific_log2_fc: float
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
            "contrasts_file": str(self.contrasts_file),
            "induction_log2_fc": self.induction_log2_fc,
            "specific_log2_fc": self.specific_log2_fc,
        }


def write_target_vs_other_conditions_contrasts(
    *,
    output_file: Path,
    design_conditions: Sequence[str],
    target_conditions: Sequence[str],
) -> Path:
    """Write directional pairwise and one-vs-mean contrast definitions.

    Conditions in ``design_conditions`` that are absent from ``target_conditions``
    remain in the design with coefficient zero (for example, an untreated control).
    """
    design_conditions = tuple(design_conditions)
    target_conditions = tuple(target_conditions)

    if len(design_conditions) < 2:
        raise ValueError("At least two design conditions are required.")
    if len(target_conditions) < 2:
        raise ValueError("At least two target conditions are required.")
    if len(set(design_conditions)) != len(design_conditions):
        raise ValueError("Design conditions must be unique.")
    if len(set(target_conditions)) != len(target_conditions):
        raise ValueError("Target conditions must be unique.")

    unknown_targets = set(target_conditions) - set(design_conditions)
    if unknown_targets:
        raise ValueError(f"Target conditions are absent from the design: {sorted(unknown_targets)}")

    rows: list[dict[str, str | float]] = []

    for target in target_conditions:
        other_conditions = tuple(condition for condition in target_conditions if condition != target)

        for other in other_conditions:
            coefficients = {condition: 0.0 for condition in design_conditions}
            coefficients[target] = 1.0
            coefficients[other] = -1.0
            rows.append({
                "contrast": f"{target}_vs_{other}",
                **coefficients,
            })

        # With one comparator this would duplicate the pairwise contrast.
        if len(other_conditions) > 1:
            coefficients = {condition: 0.0 for condition in design_conditions}
            coefficients[target] = 1.0
            other_weight = -1.0 / len(other_conditions)
            for other in other_conditions:
                coefficients[other] = other_weight
            rows.append({
                "contrast": f"{target}_vs_other_conditions_mean",
                **coefficients,
            })

    output_file.parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame(rows).write_csv(output_file)
    return output_file



def write_target_vs_control_and_other_targets_contrasts(
    *,
    output_file: Path,
    design_conditions: Sequence[str],
    target_conditions: Sequence[str],
    control_condition: str = "ctrl",
) -> Path:
    design_conditions = tuple(design_conditions)
    target_conditions = tuple(target_conditions)

    if control_condition not in design_conditions:
        raise ValueError(f"Control condition '{control_condition}' is absent from the design.")

    if control_condition in target_conditions:
        raise ValueError("Control condition must not be included in target conditions.")

    unknown_targets = set(target_conditions) - set(design_conditions)
    if unknown_targets:
        raise ValueError(f"Target conditions are absent from the design: {sorted(unknown_targets)}")

    rows: list[dict[str, str | float]] = []

    for target in target_conditions:
        # Require induction: target IFN > untreated control
        coefficients = {
            condition: 0.0 for condition in design_conditions
        }
        coefficients[target] = 1.0
        coefficients[control_condition] = -1.0

        rows.append({
            "contrast": f"{target}_vs_{control_condition}",
            **coefficients,
        })

        # Require specificity against every other IFN individually
        for other in target_conditions:
            if other == target:
                continue

            coefficients = {
                condition: 0.0 for condition in design_conditions
            }
            coefficients[target] = 1.0
            coefficients[other] = -1.0

            rows.append({
                "contrast": f"{target}_vs_{other}",
                **coefficients,
            })

    output_file.parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame(rows).write_csv(output_file)

    return output_file


def run_multi_conditional_limma_voom(limma_input: MultiConditionalLimmaVoomInput) -> None:
    r_script: Path = find_file_in_project("multi_conditional_limma_voom.R")

    run_rscript(
        rscript_exe=RSCRIPT_EXE,
        r_script=r_script,
        config=limma_input.to_config(),
    )
