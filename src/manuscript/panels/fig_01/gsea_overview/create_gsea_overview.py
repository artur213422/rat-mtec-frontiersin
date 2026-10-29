import math
import textwrap
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import polars as pl
from matplotlib.axes import Axes
from matplotlib.colors import to_rgb
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

from src.common.paths import PROJECT_ROOT
from src.manuscript.figures.layout.panel_layout import PanelLayout
from src.manuscript.shared.config.image_config import TITLE_SIZE, ManuscriptColors, DPI
from src.manuscript.shared.fig.fig import split_axis, create_figure, save_svg

MAIN_LABEL_SIZE = 17
TERM_LABEL_ALPHA = 0.23
SATURATION = 0.32

NS = {
    "svg": "http://www.w3.org/2000/svg"
}

DOWNREGULATED_FILE: Path = PROJECT_ROOT / "data/mtec_high/derived_data/standard_gsea/rrvgo/downregulated/rat_mtec_high_BP/output_reduced_terms.csv"
UPREGULATED_FILE: Path = PROJECT_ROOT / "data/mtec_high/derived_data/standard_gsea/rrvgo/upregulated/rat_mtec_high_BP/output_reduced_terms.csv"

DOWNREGULATED_SVG: Path = PROJECT_ROOT / "data/mtec_high/derived_data/standard_gsea/rrvgo/downregulated/rat_mtec_high_BP/output_treemap.svg"
UPREGULATED_SVG: Path = PROJECT_ROOT / "data/mtec_high/derived_data/standard_gsea/rrvgo/upregulated/rat_mtec_high_BP/output_treemap.svg"

def _box(element):
    return tuple(
        float(element.get(k))
        for k in ("x", "y", "width", "height")
    )


def _contains(box, x, y):
    bx, by, width, height = box

    return (
        bx - 0.02 <= x <= bx + width + 0.02
        and
        by - 0.02 <= y <= by + height + 0.02
    )


def _style(element):
    return dict(
        tuple(value.strip() for value in part.split(":", 1))
        for part in element.get("style", "").split(";")
        if ":" in part
    )


def create_rat_bp_treemap(ax: Axes, terms: pl.DataFrame, svg_file: Path) -> dict:
    if ax.has_data():
        raise ValueError("Treemap requires an empty Axes.")

    rows = terms.to_dicts()

    if not rows or len({row["go"] for row in rows}) != len(rows):
        raise ValueError("Reduced GO terms must be nonempty and unique.")

    groups = defaultdict(list)

    for row in rows:
        score = float(row["score"])

        if not math.isfinite(score) or score <= 0:
            raise ValueError(f"Invalid score for {row['go']}")

        groups[row["parentTerm"]].append(row)

    root = ET.parse(svg_file).getroot()

    leaves = []
    group_boxes = []

    for rect in root.findall(".//svg:rect", NS):

        style = _style(rect)

        stroke = style.get(
            "stroke-width",
            "",
        ).strip()

        if stroke == "2.25":
            group_boxes.append(_box(rect))

        elif stroke == "0.75" and style.get("fill", "").startswith("#"):
            leaves.append(
                (
                    _box(rect),
                    style["fill"].strip(),
                )
            )

    if len(leaves) != len(rows) or len(group_boxes) != len(groups):
        raise ValueError(
            f"SVG geometry does not match dataframe: {svg_file} has "
            f"{len(leaves)} leaves and {len(group_boxes)} groups; "
            f"expected {len(rows)} leaves and {len(groups)} groups."
        )

    # Recover group labels
    group_lines = defaultdict(list)
    term_texts = []

    for text in root.findall(".//svg:text", NS):

        style = _style(text)

        x = float(text.get("x"))
        y = float(text.get("y"))

        if style.get("font-weight", "").strip() == "bold":

            containing = [
                i
                for i, box in enumerate(group_boxes)
                if _contains(box, x, y)
            ]

            if len(containing) != 1:
                raise ValueError("Cannot associate group label.")

            group_lines[containing[0]].append(
                (
                    y,
                    text.text or "",
                )
            )

        else:
            term_texts.append(
                (
                    x,
                    y,
                    text.text or "",
                    style,
                )
            )

    names = [
        " ".join(
            line
            for _, line
            in sorted(group_lines[i])
        )
        for i in range(len(group_boxes))
    ]

    if set(names) != set(groups):
        raise ValueError("SVG labels do not match dataframe.")

    # Validate areas
    area = sum(
        width * height
        for _, _, width, height
        in group_boxes
    )

    total_score = sum(
        float(row["score"])
        for row in rows
    )

    for name, (_, _, width, height) in zip(names, group_boxes):

        expected = (
            sum(
                float(row["score"])
                for row in groups[name]
            )
            / total_score
        )

        actual = width * height / area

        if not math.isclose(
            actual,
            expected,
            rel_tol=0.003,
        ):
            raise ValueError(f"Area mismatch for {name}")

    xmin = min(
        box[0]
        for box in group_boxes
    )

    ymin = min(
        box[1]
        for box in group_boxes
    )

    xmax = max(
        box[0] + box[2]
        for box in group_boxes
    )

    ymax = max(
        box[1] + box[3]
        for box in group_boxes
    )

    ax.set_xlim(xmin, xmax)

    ax.set_ylim(ymax, ymin)

    ax.figure.canvas.draw()

    extent = ax.get_window_extent()

    source_scale = min(
        extent.width / ax.figure.dpi * 72 / (xmax - xmin),
        extent.height / ax.figure.dpi * 72 / (ymax - ymin),
    )

    # Draw leaves
    for box, color in leaves:

        rgb = to_rgb(color)

        muted = tuple(
            0.97 * (1 - SATURATION)
            + channel * SATURATION
            for channel in rgb
        )

        ax.add_patch(
            Rectangle(
                box[:2],
                box[2],
                box[3],
                facecolor=muted,
                edgecolor=ManuscriptColors.WHITE,
                linewidth=0.35 * source_scale,
            )
        )

    # Small labels
    for x, y, label, style in term_texts:

        size = float(style.get("font-size", "11px").replace("px", ""))

        ax.text(
            x,
            y,
            label,
            ha="center",
            va="baseline",
            fontsize=size * 0.8 * source_scale,
            color=ManuscriptColors.BLACK,
            alpha=TERM_LABEL_ALPHA,
            clip_on=True,
        )

    ax.set_axis_off()
    ax.set_facecolor(ManuscriptColors.WHITE)

    ax.figure.canvas.draw()

    renderer = ax.figure.canvas.get_renderer()

    # Group labels
    for name, box in zip(names, group_boxes):

        x, y, width, height = box

        ax.add_patch(
            Rectangle(
                (x, y),
                width,
                height,
                fill=False,
                edgecolor=ManuscriptColors.WHITE,
                linewidth=1.8 * source_scale,
            )
        )

        p0 = ax.transData.transform((x, y))

        p1 = ax.transData.transform((x + width, y + height))

        available_w = abs(p1[0] - p0[0]) * 0.88

        available_h = abs(p1[1] - p0[1]) * 0.82

        best = None

        for wrap_width in range(min(8, len(name)), min(len(name), 55) + 1):

            label = textwrap.fill(
                name,
                width=wrap_width,
                break_long_words=False,
                break_on_hyphens=False,
            )

            text = ax.text(
                x + width / 2,
                y + height / 2,
                label,
                ha="center",
                va="center",
                fontsize=MAIN_LABEL_SIZE,
                fontweight="bold",
            )

            extent = text.get_window_extent(renderer)

            size = (
                MAIN_LABEL_SIZE
                *
                min(
                    1,
                    available_w / extent.width,
                    available_h / extent.height,
                )
            )

            text.remove()

            if best is None or size > best[0]:
                best = (size, label)

        size, label = best

        ax.text(
            x + width / 2,
            y + height / 2,
            label,
            ha="center",
            va="center",
            fontsize=size,
            fontweight="bold",
            color="#25313b",
            bbox={
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.60,
                "pad": 2 * source_scale,
            },
        )


    return {
        "terms": len(rows),
        "groups": len(groups),
        "leaf_area": area,
    }


def create_gsea_overview(layout: PanelLayout, output_file: Path) -> Path:
    fig: Figure = create_figure(width_mm=layout.width, height_mm=layout.height, dpi=DPI)

    ax = fig.add_axes((0, 0, 1, 1))

    ax_down, ax_up = split_axis(ax, rows=1,cols=2,)

    down_terms: pl.DataFrame = pl.read_csv(DOWNREGULATED_FILE)
    up_terms: pl.DataFrame = pl.read_csv(UPREGULATED_FILE)

    create_rat_bp_treemap(ax_down, down_terms,DOWNREGULATED_SVG)
    create_rat_bp_treemap(ax_up, up_terms, UPREGULATED_SVG)

    ax_down.set_title("Downregulated", fontsize=TITLE_SIZE, color=ManuscriptColors.BLACK)
    ax_up.set_title("Upregulated", fontsize=TITLE_SIZE, color=ManuscriptColors.BLACK)

    save_svg(fig, output_file)

    return output_file