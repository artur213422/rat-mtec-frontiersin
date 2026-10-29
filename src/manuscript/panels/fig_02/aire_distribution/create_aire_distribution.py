from __future__ import annotations

from importlib import import_module
from pathlib import Path
import subprocess
import sys
from matplotlib.figure import Figure

from src.common.paths import PROJECT_ROOT
from src.manuscript.figures.layout.panel_layout import PanelLayout
from src.manuscript.shared.config.image_config import DPI
from src.manuscript.shared.fig.fig import create_figure, save_svg

PACKAGE = "src.manuscript.panels.fig_02.aire_distribution"

STEPS = (
    "01_pick_aire_top_genes",
    "02_count_aire_neighbors",
    "03_random_nearest_neighbor_permutations",
)

visualize = import_module(f"{PACKAGE}.steps.04_visualize").visualize

def create_aire_distribution(layout: PanelLayout, output_file: Path, recalculate: bool = False) -> Path:

    if recalculate:
        for index, step in enumerate(STEPS, start=1):

            print(f"\n[{index}/{len(STEPS)}] Running {step}", flush=True,)

            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    f"{PACKAGE}.steps.{step}",
                ],
                cwd=PROJECT_ROOT,
                check=True,
            )

    fig: Figure = create_figure(width_mm=layout.width, height_mm=layout.height, dpi=DPI,)

    ax = fig.add_axes((0, 0, 1, 1))

    visualize(ax)

    save_svg(fig, output_file)

    return output_file