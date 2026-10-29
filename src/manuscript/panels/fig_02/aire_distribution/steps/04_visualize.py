from matplotlib.axes import Axes
from matplotlib.lines import Line2D

import numpy as np
import numpy.typing as npt
import polars as pl
from scipy.stats import gaussian_kde

from src.manuscript.panels.fig_02.aire_distribution.confg import PERMUTATIONS_FILE, NEIGHBOR_COUNTS_FILE, \
    NEIGHBOR_DISTANCE_KB
from src.manuscript.shared.config.image_config import (
    FONT_SIZE,
    SCALE_SIZE,
    X_AXIS_SIZE,
    ManuscriptColors,
)


DATASET_GROUPS: dict[str, str] = {
    "rat_mtec_high": "Rat",
    "mouse_gse85_mtec_high": "B6",
    "mouse_gse8564_b6_mtec_high": "B6",
    "mouse_gse2585_mtec_high": "B6",
    "mouse_gse33878_mtec_high": "NOD",
    "mouse_gse8564_nod_mtec_high": "NOD",
    "mouse_gse8564_balbc_mtec_high": "BALB/c",
}

GROUP_COLORS: dict[str, str] = {
    "Rat": ManuscriptColors.RED,
    "B6": ManuscriptColors.BLACK,
    "NOD": ManuscriptColors.GREEN,
    "BALB/c": ManuscriptColors.BLUE,
}

RANDOM_DISTRIBUTION_COLOR: str = "gray"

EFFECT_STYLES: dict[str, str] = {
    "activated": "-",
    "repressed": "--",
}


def visualize(ax: Axes) -> None:
    permutation_df: pl.DataFrame = pl.read_csv(PERMUTATIONS_FILE)
    observed_df: pl.DataFrame = pl.read_csv(NEIGHBOR_COUNTS_FILE)

    # ------------------------------------------------------------------
    # Plot distributions and observed values
    # ------------------------------------------------------------------

    dataset_ids: list[str] = (
        observed_df
        .get_column("dataset_id")
        .unique(maintain_order=True)
        .to_list()
    )

    for dataset_id in dataset_ids:
        group: str = DATASET_GROUPS[dataset_id]
        color: str = GROUP_COLORS[group]

        dataset_perm_df: pl.DataFrame = permutation_df.filter(
            pl.col("dataset_id") == dataset_id
        )

        dataset_observed_df: pl.DataFrame = observed_df.filter(
            pl.col("dataset_id") == dataset_id
        )

        plotted_sizes: set[int] = set()

        for effect, linestyle in EFFECT_STYLES.items():
            effect_df: pl.DataFrame = dataset_observed_df.filter(
                pl.col("aire_effect") == effect
            )

            if effect_df.is_empty():
                continue

            if effect_df.height != 1:
                raise ValueError(
                    f"Duplicate observed counts for {dataset_id}/{effect}."
                )

            effect_perm_df: pl.DataFrame = dataset_perm_df.filter(
                pl.col("aire_effect") == effect
            )

            n_genes: int = int(effect_df["n_genes"][0])

            if effect_perm_df.is_empty():
                raise ValueError(
                    f"Missing permutations for {dataset_id}/{effect}."
                )

            if effect_perm_df["n_genes"].unique().to_list() != [n_genes]:
                raise ValueError(
                    f"Permutation gene counts disagree for {dataset_id}/{effect}."
                )

            # ----------------------------------------------------------
            # Random distribution
            # ----------------------------------------------------------

            if n_genes not in plotted_sizes:
                values: npt.NDArray[np.float64] = (
                    effect_perm_df
                    .get_column("n_with_neighbor")
                    .cast(pl.Float64)
                    .to_numpy()
                )

                if values.size < 2 or not np.isfinite(values).all():
                    raise ValueError(
                        f"Invalid permutation values for {dataset_id}/{effect}."
                    )

                if np.ptp(values) == 0:
                    ax.axvline(
                        values[0],
                        color=RANDOM_DISTRIBUTION_COLOR,
                        alpha=0.4,
                        linewidth=2,
                    )
                else:
                    grid: npt.NDArray[np.float64] = np.linspace(
                        max(0, float(values.min()) - 5),
                        float(values.max()) + 5,
                        500,
                    )

                    ax.plot(
                        grid,
                        gaussian_kde(values)(grid),
                        color=RANDOM_DISTRIBUTION_COLOR,
                        alpha=0.4,
                        linewidth=2,
                        linestyle=linestyle,
                    )

                plotted_sizes.add(n_genes)

            # ----------------------------------------------------------
            # Observed values
            # ----------------------------------------------------------

            ax.axvline(
                int(effect_df["n_with_neighbor"][0]),
                color=color,
                linestyle=linestyle,
                linewidth=2.5,
            )

    # ------------------------------------------------------------------
    # Axes
    # ------------------------------------------------------------------

    ax.set_xlabel(
        f"Genes with another selected gene within {NEIGHBOR_DISTANCE_KB} kb",
        fontsize=X_AXIS_SIZE,
        color=ManuscriptColors.BLACK,
    )

    ax.set_ylabel(
        "Density",
        fontsize=FONT_SIZE,
        color=ManuscriptColors.BLACK,
    )

    ax.tick_params(
        axis="both",
        labelsize=SCALE_SIZE,
        colors=ManuscriptColors.BLACK,
    )

    ax.set_facecolor(ManuscriptColors.WHITE)

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    # ------------------------------------------------------------------
    # Legend
    # ------------------------------------------------------------------

    legend_handles: list[Line2D] = [
        Line2D(
            [],
            [],
            color=color,
            linewidth=2.5,
            label=group,
        )
        for group, color in GROUP_COLORS.items()
        if any(DATASET_GROUPS[dataset_id] == group for dataset_id in dataset_ids)
    ]

    legend_handles.extend([
        Line2D(
            [],
            [],
            color=ManuscriptColors.BLACK,
            linestyle="-",
            linewidth=1.5,
            label="Activated",
        ),
        Line2D(
            [],
            [],
            color=ManuscriptColors.BLACK,
            linestyle="--",
            linewidth=1.5,
            label="Repressed",
        ),
        Line2D(
            [],
            [],
            color=RANDOM_DISTRIBUTION_COLOR,
            linewidth=2,
            alpha=0.4,
            label="Random distribution",
        ),
    ])

    ax.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        frameon=False,
        fontsize=FONT_SIZE,
        labelcolor=ManuscriptColors.BLACK,
    )