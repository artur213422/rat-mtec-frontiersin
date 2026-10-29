"""Join model evidence before annotation or selection; never union contrasts."""

from typing import Sequence

import polars as pl


def join_induction_and_specificity(results: pl.DataFrame, target_conditions: Sequence[str]) -> pl.DataFrame:
    target_conditions = tuple(target_conditions)
    joined_targets: list[pl.DataFrame] = []

    for target in target_conditions:
        required_contrasts = [
            ("induction", f"{target}_vs_ctrl"),
            *[
                (f"specific_vs_{other}", f"{target}_vs_{other}")
                for other in target_conditions
                if other != target
            ],
        ]

        evidence_frames: list[pl.DataFrame] = []

        for prefix, contrast in required_contrasts:
            frame = results.filter(pl.col("contrast") == contrast)

            if frame.is_empty():
                raise ValueError(f"Missing contrast {contrast}; rerun the limma step.")

            if frame["gene_symbol"].null_count() or frame["gene_symbol"].n_unique() != frame.height:
                raise ValueError(f"Non-unique or missing gene symbols in {contrast}")

            frame = frame.rename({
                column: f"{prefix}_{column}"
                for column in frame.columns
                if column != "gene_symbol"
            })

            evidence_frames.append(frame)

        joined = evidence_frames[0]

        for frame in evidence_frames[1:]:
            if set(joined["gene_symbol"]) != set(frame["gene_symbol"]):
                raise ValueError(f"Different gene universes for {target}")

            joined = joined.join(
                frame,
                on="gene_symbol",
                validate="1:1",
            )

        joined_targets.append(
            joined.with_columns(
                pl.lit(target).alias("target_condition")
            )
        )

    return pl.concat(
        joined_targets,
        how="diagonal",
    ).sort(
        ["target_condition", "gene_symbol"]
    )


def pick_induced_specific(
    candidates: pl.DataFrame,
    *,
    target_conditions: Sequence[str],
    induction_log2_fc: float = 1.0,
    specific_log2_fc: float = 0.5,
    fdr: float = 0.05,
) -> pl.DataFrame:
    target_conditions = tuple(target_conditions)
    selected_targets: list[pl.DataFrame] = []

    for target in target_conditions:
        target_df = candidates.filter(
            pl.col("target_condition") == target
        )

        specificity_rules: list[pl.Expr] = []

        for other in target_conditions:
            if other == target:
                continue

            specificity_rules.extend([
                (
                    pl.col(
                        f"specific_vs_{other}_log2_fold_change"
                    )
                    >= specific_log2_fc
                ),
                (
                    pl.col(
                        f"specific_vs_{other}_p_adj"
                    )
                    <= fdr
                ),
            ])

        selected = target_df.filter(
            # Target IFN is induced versus control
            (
                pl.col("induction_log2_fold_change")
                >= induction_log2_fc
            )
            & (pl.col("induction_p_adj") <= fdr)

            # Target IFN exceeds every competing IFN
            & pl.all_horizontal(specificity_rules)
        )

        selected_targets.append(selected)

    return pl.concat(selected_targets).sort(
        ["target_condition", "induction_log2_fold_change"],
        descending=[False, True],
    )