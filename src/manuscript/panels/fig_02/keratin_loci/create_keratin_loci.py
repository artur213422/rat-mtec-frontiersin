from pathlib import Path

import numpy as np
import polars as pl
from matplotlib.axes import Axes
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.figure import Figure

from src.common.dataset_catalog.api.get_datasets import get_rat_tec_high, get_mouse_tec_high
from src.manuscript.figures.layout.panel_layout import PanelLayout
from src.manuscript.shared.config.image_config import DPI, FONT_SIZE, SCALE_SIZE, TITLE_SIZE, X_AXIS_SIZE, ManuscriptColors
from src.manuscript.shared.fig.fig import create_figure, save_svg

# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------
mouse_dataset_id: str = "mouse_gse8564_b6_mtec_high"

mouse_input = {
    "dataset": next(dataset for dataset in get_mouse_tec_high() if dataset.dataset_id == mouse_dataset_id),
    "chromosomes": ["chr11", "chr15"],
}

rat_input = {
    "dataset": get_rat_tec_high(),
    "chromosomes": ["chr10", "chr7"]
}

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
def _draw_loci(keratins: list[pl.DataFrame], ax: Axes) -> None:
    if not keratins or any(data.is_empty() for data in keratins):
        raise ValueError("Each keratin locus must contain genes")

    loci = []

    for data in keratins:
        if data["chromosome"].null_count() or data["chromosome"].n_unique() != 1:
            raise ValueError("Each keratin locus must belong to one chromosome")

        data = data.sort("start")

        x = data["start"].to_numpy().astype(float) / 1_000_000
        y = data["log2_fold_change"].to_numpy().astype(float)

        if not np.isfinite(x).all() or not np.isfinite(y).all():
            raise ValueError("Keratin genes must have finite starts and fold changes")
        loci.append((data, x, y))

    limit = max(1.0, *(float(np.abs(y).max()) for _, _, y in loci))
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)

    cmap = LinearSegmentedColormap.from_list(
        "keratin_fc",
        [ManuscriptColors.BLUE, ManuscriptColors.WHITE, ManuscriptColors.RED],
    )

    figure = ax.figure
    grid = ax.get_subplotspec().subgridspec(1, len(loci))
    ax.remove()
    shared_axis = None

    for index, (data, x, y) in enumerate(loci):
        locus_ax = figure.add_subplot(grid[0, index], sharey=shared_axis)

        if shared_axis is None:
            shared_axis = locus_ax

        p_values = data["p_val"].to_numpy().astype(float)
        significant = np.isfinite(p_values) & (p_values < 0.05)
        locus_ax.plot(x, y, color=ManuscriptColors.BLACK, alpha=0.3, linewidth=0.8)
        locus_ax.scatter(
            x[significant], y[significant], c=y[significant],
            cmap=cmap, norm=norm, s=20, zorder=3,
        )

        locus_ax.scatter(
            x[~significant], y[~significant], color=ManuscriptColors.BLACK,
            alpha=0.4, s=20, zorder=3,
        )

        locus_ax.axhline(0, color=ManuscriptColors.BLACK, alpha=0.4, linewidth=0.5)

        for threshold in (-1, 1):
            locus_ax.axhline(
                threshold, color=ManuscriptColors.RED, linestyle="--",
                linewidth=0.6, alpha=0.6,
            )

        for symbol, position, fold_change in zip(data["gene_symbol"], x, y):
            locus_ax.annotate(
                symbol, (position, fold_change), xytext=(2, 5),
                textcoords="offset points", rotation=90, ha="left", va="bottom",
                fontsize=FONT_SIZE, color=ManuscriptColors.BLACK,
            )

        padding = max(float(np.ptp(x)) * 0.08, 0.005)
        locus_ax.set_xlim(float(x.min()) - padding, float(x.max()) + padding)
        locus_ax.set_ylim(-limit * 1.3, limit * 1.3)
        locus_ax.set_title(data["chromosome"][0], fontsize=TITLE_SIZE)
        locus_ax.set_xlabel("Gene start (Mb)", fontsize=FONT_SIZE)

        if index == 0:
            locus_ax.set_ylabel(r"log$_2$ fold change", fontsize=X_AXIS_SIZE)

        locus_ax.tick_params(
            labelsize=SCALE_SIZE, colors=ManuscriptColors.BLACK,
            labelleft=index == 0,
        )

        locus_ax.ticklabel_format(axis="x", style="plain", useOffset=False)

        for side in ("top", "right"):
            locus_ax.spines[side].set_visible(False)

        for side in ("left", "bottom"):
            locus_ax.spines[side].set_color(ManuscriptColors.BLACK)

def _prepare_keratin_loci(tec_df: pl.DataFrame, chromosomes: list[str]) -> list[pl.DataFrame]:
    keratins = tec_df.filter(pl.col("gene_symbol").str.contains(r"(?i)^krt\d+[a-z]*$"))

    return [
        keratins.filter(pl.col("chromosome") == chromosome).sort("start")
        for chromosome in chromosomes
    ]

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def create_keratin_loci(layout: PanelLayout, output_file: Path) -> Path:
    fig: Figure = create_figure(width_mm=layout.width, height_mm=layout.height, dpi=DPI)

    fig.set_layout_engine("constrained")
    axes = fig.subplots(2, 1)

    ## Rat
    rat_keratins: list[pl.DataFrame] =_prepare_keratin_loci(
        tec_df=rat_input["dataset"].get_standardized_df(),
        chromosomes=rat_input["chromosomes"]
    )

    _draw_loci(rat_keratins, axes[0])

    ## Mouse
    mouse_keratins: list[pl.DataFrame] = _prepare_keratin_loci(
        tec_df=mouse_input["dataset"].get_standardized_df(),
        chromosomes=mouse_input["chromosomes"]
    )

    _draw_loci(mouse_keratins, axes[1])

    save_svg(fig, output_file)

    return output_file
