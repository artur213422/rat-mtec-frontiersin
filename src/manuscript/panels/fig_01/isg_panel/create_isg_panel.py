import math
from pathlib import Path

import polars as pl
from matplotlib.axes import Axes
from matplotlib.lines import Line2D

from src.manuscript.figures.layout.panel_layout import PanelLayout
from src.manuscript.panels.fig_01.isg_panel._read_files import (
    IFN_VS_CTRL_FILE,
    INDUCED_ISG_FILE,
)
from src.manuscript.shared.config.dataset_labels import DATASET_LABELS
from src.manuscript.shared.config.image_config import (
    FONT_SIZE,
    SCALE_SIZE,
    TITLE_SIZE,
    ManuscriptColors,
    DPI,
)
from src.manuscript.shared.fig.fig import create_figure, save_svg


FDR_THRESHOLD = 0.05
NES_RANGE = (-3.0, 1.75)
MARKER_SIZE = 25

SERIES = (
    ("ISG response to IFN", ManuscriptColors.BLUE, "o"),
    ("IFN-specific ISG induction", ManuscriptColors.RED, "s"),
)

def _read_results(path: Path, term: str) -> pl.DataFrame:
    return (
        pl.read_csv(path, null_values=["NA", ""])
        .filter(
            (pl.col("term") == term)
            & pl.col("dataset_id").is_in(list(DATASET_LABELS))
        )
        .select(
            pl.col("dataset_id").alias("dataset"),
            pl.col("nes").cast(pl.Float64),
            pl.col("fdr_q_value").cast(pl.Float64).alias("fdr"),
        )
    )


def _get_result(
    data: pl.DataFrame,
    dataset_id: str,
    series_label: str,
) -> tuple[float, float]:
    row = data.filter(pl.col("dataset") == dataset_id)

    if row.height != 1:
        raise ValueError(
            f"{series_label}: expected one result for {dataset_id}, "
            f"found {row.height}"
        )

    nes = row["nes"][0]
    fdr = row["fdr"][0]

    if nes is None or not math.isfinite(nes):
        raise ValueError(
            f"{series_label}: missing or non-finite NES for {dataset_id}"
        )

    return float(nes), float(fdr) if fdr is not None else math.nan


def _create_dumbbell_graph(
    ax: Axes,
    ifn_vs_ctrl: pl.DataFrame,
    induced: pl.DataFrame,
) -> None:
    for y, dataset_id in enumerate(DATASET_LABELS):
        results = [
            _get_result(data, dataset_id, label)
            for data, (label, _, _) in zip(
                (ifn_vs_ctrl, induced), SERIES
            )
        ]

        ax.plot(
            [result[0] for result in results],
            [y, y],
            color=ManuscriptColors.BLACK,
            alpha=0.35,
            linewidth=1.2,
            zorder=2,
        )

        for (nes, fdr), (_, color, marker) in zip(results, SERIES):
            significant = math.isfinite(fdr) and fdr < FDR_THRESHOLD

            ax.scatter(
                nes,
                y,
                s=MARKER_SIZE,
                marker=marker,
                facecolors=color if significant else ManuscriptColors.WHITE,
                edgecolors=color,
                linewidths=1,
                zorder=3,
            )

    ax.axvline(
        0,
        color=ManuscriptColors.BLACK,
        alpha=0.3,
        linewidth=0.8,
        zorder=0,
    )

    ax.set_xlim(*NES_RANGE)
    ax.set_ylim(len(DATASET_LABELS) - 0.5, -0.5)
    ax.set_yticks(
        range(len(DATASET_LABELS)),
        list(DATASET_LABELS.values()),
        fontsize=FONT_SIZE,
        color=ManuscriptColors.BLACK,
    )
    ax.set_xticks([-3, -2, -1, 0, 1])
    ax.set_xlabel("NES", fontsize=FONT_SIZE)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(
        axis="x",
        labelsize=SCALE_SIZE,
        colors=ManuscriptColors.BLACK,
    )
    ax.grid(
        axis="x",
        color=ManuscriptColors.BLACK,
        alpha=0.08,
        linewidth=0.5,
    )
    ax.set_axisbelow(True)

    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)

    ax.spines["bottom"].set_color(ManuscriptColors.BLACK)


def create_isg_panel(layout: PanelLayout, output_file: Path) -> Path:
    fig = create_figure(
        width_mm=layout.width,
        height_mm=layout.height,
        dpi=DPI,
    )
    fig.set_layout_engine("constrained")
    axes = fig.subplots(2, 1)

    for ax, term, title in zip(
        axes,
        ("ifna", "ifng"),
        ("Type I ISGs", "Type II ISGs"),
    ):
        _create_dumbbell_graph(
            ax,
            _read_results(IFN_VS_CTRL_FILE, term),
            _read_results(INDUCED_ISG_FILE, term),
        )
        ax.set_title(
            title,
            fontsize=TITLE_SIZE,
            color=ManuscriptColors.BLACK,
        )

    handles = [
        Line2D(
            [],
            [],
            linestyle="none",
            marker=marker,
            markerfacecolor=color,
            markeredgecolor=color,
            markersize=math.sqrt(MARKER_SIZE),
            label=label,
        )
        for label, color, marker in SERIES
    ]

    axes[-1].legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.35),
        frameon=False,
        fontsize=FONT_SIZE,
        title="Filled: FDR < 0.05; open: otherwise",
        title_fontsize=SCALE_SIZE,
    )

    save_svg(fig, output_file)
    return output_file