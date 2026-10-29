"""Compare rat mTEC low/high GSEA enrichment of cell-type signatures."""

from pathlib import Path
from matplotlib import pyplot as plt
from src.manuscript.figures.layout.panel_layout import PanelLayout
from src.manuscript.shared.fig.fig import create_figure, save_svg, add_panel_label
from src.manuscript.shared.config.image_config import DPI

import csv
import math

from matplotlib.axes import Axes
from matplotlib.lines import Line2D

from src.common.paths import PROJECT_ROOT
from src.manuscript.shared.config.image_config import (
    FONT_SIZE, SCALE_SIZE, TITLE_SIZE, X_AXIS_SIZE, ManuscriptColors,
)

SOURCES = [
    ("Michelson et al., 2022", "michelson_2022_mimetic_cells", "gsea_analysis.csv"),
    ("Baran-Gale et al., 2020", "hollander_thymus_cell_types_gsea", "cell_type_gsea.csv"),
]
POPULATIONS = {
    "mtec_low": ("mTEC low", ManuscriptColors.BLUE, "o"),
    "mtec_high": ("mTEC high", ManuscriptColors.RED, "s"),
}
FDR_THRESHOLD = 0.05
NES_LIMIT = 2.5
MARKER_SIZE = 22
SECTION_GAP = 0.8

IGNORE_POPULATIONS: dict[str, list[str]] = {
    "Michelson et al., 2022": [],
    "Baran-Gale et al., 2020": ["Mature cTEC", "Perinatal cTEC"],
}

def _load_rat_results(folder: str, filename: str, population: str) -> dict:
    path = PROJECT_ROOT / "data" / population / "derived_data" / folder / filename
    results = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"dataset_id", "term", "nes", "fdr_q_value"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
        for row in reader:
            if row["dataset_id"] != f"rat_{population}":
                continue
            term = row["term"]
            if term in results:
                raise ValueError(f"Duplicate rat signature {term!r} in {path}")
            values = []
            for column in ("nes", "fdr_q_value"):
                raw = row[column].strip()
                values.append(float(raw) if raw and raw.upper() != "NA" else math.nan)
            results[term] = tuple(values)
    if not results:
        raise ValueError(f"No rat_{population} results found in {path}")
    return results


def _draw_panel(ax: Axes) -> None:
    """Draw paired NES values; filled markers indicate FDR < 0.05."""
    if ax.has_data():
        raise ValueError("Mimetic dumbbell plot requires an empty Axes.")

    # Read and validate everything before drawing, keeping identical names from
    # different signature collections separate.
    sections = []
    for title, folder, filename in SOURCES:
        results = {
            population: _load_rat_results(folder, filename, population)
            for population in POPULATIONS
        }

        ignored = set(IGNORE_POPULATIONS.get(title, []))
        results = {
            population: {
                term: result
                for term, result in values.items()
                if term not in ignored
            }
            for population, values in results.items()
        }

        terms = sorted(
            set().union(*(set(values) for values in results.values()))
        )
        if not terms:
            continue

        for population, values in results.items():
            missing = set(terms) - set(values)
            if missing:
                raise ValueError(
                    f"Missing {population} signatures in {title}: "
                    f"{sorted(missing)}"
                )

        for values in results.values():
            for term, (nes, _) in values.items():
                if math.isfinite(nes) and abs(nes) > NES_LIMIT:
                    raise ValueError(
                        f"NES for {title}/{term} exceeds "
                        f"NES_LIMIT={NES_LIMIT}."
                    )

        sections.append((title, terms, results))

    positions, labels = [], []
    cursor = 0.0
    ax.axvline(0, color=ManuscriptColors.BLACK, alpha=0.3, linewidth=0.7, zorder=0)
    for title, terms, results in sections:
        ax.text(0, cursor, title, transform=ax.get_yaxis_transform(),
                ha="left", va="center", fontsize=TITLE_SIZE,
                fontweight="bold", color=ManuscriptColors.BLACK)
        cursor += 1
        section_positions = []
        for term in terms:
            positions.append(cursor)
            section_positions.append(cursor)
            labels.append(term)
            low = results["mtec_low"][term][0]
            high = results["mtec_high"][term][0]
            if math.isfinite(low) and math.isfinite(high):
                ax.plot([low, high], [cursor, cursor],
                        color=ManuscriptColors.BLACK, alpha=0.25,
                        linewidth=1, zorder=1)
            cursor += 1

        # Draw squares first so nearby circle/square pairs stay distinguishable.
        for population in reversed(POPULATIONS):
            _, color, marker = POPULATIONS[population]
            for significant in (False, True):
                points = [
                    (results[population][term][0], y)
                    for term, y in zip(terms, section_positions)
                    if math.isfinite(results[population][term][0])
                    and (math.isfinite(results[population][term][1])
                         and results[population][term][1] < FDR_THRESHOLD) == significant
                ]
                if points:
                    x, y = zip(*points)
                    ax.scatter(x, y, s=MARKER_SIZE, marker=marker,
                               facecolors=color if significant else ManuscriptColors.WHITE,
                               edgecolors=color, linewidths=0.9, zorder=3)
        cursor += SECTION_GAP

    ax.set_xlim(-NES_LIMIT, NES_LIMIT)
    ax.set_ylim(cursor - SECTION_GAP - 0.4, -0.7)
    ax.set_yticks(positions, labels, fontsize=FONT_SIZE, color=ManuscriptColors.BLACK)
    ax.set_xticks([-2, -1, 0, 1, 2])
    ax.set_xlabel("NES (KO versus control)", fontsize=X_AXIS_SIZE,
                  color=ManuscriptColors.BLACK)
    ax.tick_params(axis="x", labelsize=SCALE_SIZE, colors=ManuscriptColors.BLACK)
    ax.tick_params(axis="y", length=0)
    ax.set_axisbelow(True)
    ax.grid(axis="x", color=ManuscriptColors.BLACK, alpha=0.08, linewidth=0.5)
    ax.set_facecolor(ManuscriptColors.WHITE)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(ManuscriptColors.BLACK)

    handles = [
        Line2D([], [], linestyle="none", marker=marker, color=color,
               markerfacecolor=color, markersize=4.5, label=label)
        for label, color, marker in POPULATIONS.values()
    ]
    handles.extend([
        Line2D([], [], linestyle="none", marker="o", color=ManuscriptColors.BLACK,
               markerfacecolor=ManuscriptColors.BLACK, markersize=4.5,
               label=f"FDR < {FDR_THRESHOLD:g}"),
        Line2D([], [], linestyle="none", marker="o", color=ManuscriptColors.BLACK,
               markerfacecolor=ManuscriptColors.WHITE, markersize=4.5,
               label="Not significant / FDR missing"),
    ])
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.12),
              ncol=2, frameon=False, fontsize=FONT_SIZE,
              labelcolor=ManuscriptColors.BLACK, columnspacing=1, handletextpad=0.4)



def create_mimetic_panel(layout: PanelLayout, output_file: Path) -> Path:
    """Build a standalone panel SVG for manuscript composition."""
    fig = create_figure(width_mm=layout.width, height_mm=layout.height, dpi=DPI)
    try:
        fig.set_layout_engine("constrained")
        ax = fig.subplots()
        _draw_panel(ax)
        add_panel_label(ax, layout.name)
        return save_svg(fig, output_file)
    finally:
        plt.close(fig)
