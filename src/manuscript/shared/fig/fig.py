from matplotlib import pyplot as plt
from pathlib import Path
from matplotlib.axes import Axes
from src.manuscript.shared.config.image_config import PANEL_LABEL_SIZE


plt.rcParams["font.family"] = "Arial"

def save_svg(fig: plt.Figure, output_file: Path, *,tight: bool = True) -> Path:
    try:
        output_file.parent.mkdir(parents=True, exist_ok=True)

        with plt.rc_context({"svg.fonttype": "none"}):
            fig.savefig(
                output_file,
                format="svg",
                bbox_inches="tight" if tight else None,
                pad_inches=0.04,
            )
    finally:
        plt.close(fig)
    return output_file

def create_figure(width_mm: float, height_mm: float, dpi: int = 300) -> plt.Figure:
    cm_to_inch = 1 / 2.54

    fig = plt.figure(figsize=(width_mm / 10 * cm_to_inch, height_mm / 10 * cm_to_inch), dpi=dpi)
    return fig


def add_panel_label(ax: Axes, label: str) -> None:
    ax.text(
        0.0,
        1.0,
        label,
        transform=ax.transAxes,
        fontsize=PANEL_LABEL_SIZE,
        fontweight="bold",
        va="bottom",
        ha="left",
        clip_on=False,
    )

def split_axis(ax: Axes, rows: int = 1, cols: int = 1, gap: float = 0.02) -> list[Axes]:

    axes = []

    width = (1 - gap * (cols - 1)) / cols
    height = (1 - gap * (rows - 1)) / rows

    for row in range(rows):
        for col in range(cols):

            x = col * (width + gap)
            y = 1 - (row + 1) * height - row * gap

            child = ax.inset_axes(
                (x, y, width, height)
            )

            axes.append(child)

    ax.set_axis_off()

    return axes
