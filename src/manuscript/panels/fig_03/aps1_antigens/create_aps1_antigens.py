
from pathlib import Path
import warnings

import polars as pl
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.lines import Line2D

from src.common.paths import PROJECT_ROOT
from src.manuscript.figures.layout.panel_layout import PanelLayout
from src.manuscript.shared.config.image_config import DPI, FONT_SIZE, SCALE_SIZE, X_AXIS_SIZE, ManuscriptColors
from src.manuscript.shared.fig.fig import add_panel_label,  create_figure, save_svg

# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------

ANTIGENS_FILE = PROJECT_ROOT / "data/mtec_high/derived_data/tsa_expression/aps1_known_antigens/standardized/aps1_antigens.csv"
EXPRESSION_FILE = PROJECT_ROOT / "data/mtec_high/standardized/rat/rat_mtec_high.csv"

P_VAL_THRESHOLD = 0.05
MARKER_SIZE = 28

TISSUE_GROUPS: dict[str, tuple[str, ...]] = {
    "Pancreas": ("INS", "GAD2", "PTPRN", "PNLIP", "PDX1"),
    "Endocrine organs": ("TPO", "TG", "NLRP5", "CYP11A1", "CYP17A1", "CYP21A2", "GPHB5",),
    "Gut": ("HDC", "TPH1", "RFX6", "NKX6-3", "GIP", "CDHR5"),
    "Liver": ("CYP1A2", "CYP2A6"),
    "Skin": ("SOX9", "SOX10", "TH"),
    "Neural / neuroendocrine": ("DDC", "ECE2", "ASMT", "PDYN"),
    "Lung / airway": ("BPIFB1", "KCNRG", "GAS2L2", "CROCC2"),
    "Reproductive organs": ("MAGEB2", "PDILT", "TGM4", "TSGA10", "KHDC3L", "SPAG17", "CT45A10", "ADGRG2",),
    "Dental enamel": ("ACP4", "LAMB3", "FAM20A", "ENAM", "AMELX"),
    "Kidney": ("AQP2",),
}

TISSUE_COLORS: dict[str, str] = dict(zip(
    TISSUE_GROUPS,
    ("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",),
))

SHOW_LABELS: set[str] = {
    "Pnlip", "Cyp11a1", "Hdc", "Tgm4", "Tph1", "Spag17","Tg", "Pdx1", "Cyp17a1", "Th", "Gip", "Ins2", "Bpifb1", "Sox9", "Lamb3", "Cdhr5", "Adgrg2",
}


# ------------------------------------------------------------------------------
# Prepare Data
# ------------------------------------------------------------------------------

def _load_antigens() -> pl.DataFrame:
    """Join APS-1 antigen orthologs with measured rat gene expression."""

    group_by_antigen: dict[str, str] = {}

    for group, symbols in TISSUE_GROUPS.items():
        for symbol in symbols:
            if symbol in group_by_antigen:
                raise ValueError(f"Antigen {symbol} occurs in multiple tissue groups.")
            group_by_antigen[symbol] = group

    antigens = (
        pl.read_csv(ANTIGENS_FILE, schema_overrides={"ortholog_id": pl.String})
        .filter(pl.col("dataset_id") == "rat_mtec_high")
        .select("human_gene_symbol", "ortholog_id", "tissue_related")
        .drop_nulls("ortholog_id")
        .unique(maintain_order=True)
    )

    unknown = set(antigens["human_gene_symbol"]) - set(group_by_antigen)
    if unknown:
        raise ValueError(f"Add tissue groups for antigens: {sorted(unknown)}")

    duplicates = (
        antigens.group_by("ortholog_id")
        .len()
        .filter(pl.col("len") > 1)
    )
    if not duplicates.is_empty():
        raise ValueError(f"Conflicting antigen annotations: {duplicates.to_dicts()}")

    expression = (
        pl.read_csv(EXPRESSION_FILE, schema_overrides={"ortholog_id": pl.String},)
        .filter(pl.col("ortholog_id").is_in(antigens["ortholog_id"].to_list()))
    )

    # Preserve distinct measured rat paralogs within an ortholog group.
    duplicates = (
        expression.group_by(["ortholog_id", "gene_symbol"])
        .len()
        .filter(pl.col("len") > 1)
    )
    if not duplicates.is_empty():
        raise ValueError(f"Duplicate measured rat genes: {duplicates.to_dicts()}")

    expression = expression.select(
        "ortholog_id",
        "gene_symbol",
        "p_val",
        pl.mean_horizontal("ctrl1", "ctrl2", "ctrl3").alias("ctrl_expression"),
        (pl.col("log2_fold_change")).alias("log2_fc"),
    )

    data = (
        antigens.join(expression, on="ortholog_id", how="inner")
        .with_columns(
            pl.col("human_gene_symbol")
            .replace_strict(group_by_antigen)
            .alias("tissue_group")
        )
        .filter(pl.col("ctrl_expression").is_finite() & pl.col("log2_fc").is_finite())
    )

    if data.is_empty():
        raise ValueError("No measured rat APS-1 antigens available.")

    return data


# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------

def _draw_panel(ax: Axes) -> list[Line2D]:
    if ax.has_data():
        raise ValueError("APS-1 scatter plot requires an empty Axes.")

    data = _load_antigens()

    ax.axvline(
        0,
        color=ManuscriptColors.BLACK,
        alpha=0.3,
        linewidth=0.8,
        zorder=0,
    )

    handles: list[Line2D] = []

    significant = (pl.col("p_val").is_finite() & (pl.col("p_val") < P_VAL_THRESHOLD)).fill_null(False)

    for group, color in TISSUE_COLORS.items():
        group_data = data.filter(pl.col("tissue_group") == group)

        if group_data.is_empty():
            continue

        for filled in (False, True):
            points = group_data.filter(significant if filled else ~significant)

            if points.is_empty():
                continue

            ax.scatter(
                points["log2_fc"].to_numpy(),
                points["ctrl_expression"].to_numpy(),
                s=MARKER_SIZE,
                edgecolors=color,
                facecolors=color if filled else ManuscriptColors.WHITE,
                linewidths=1,
                zorder=2,
            )

        handles.append(
            Line2D(
                [], [],
                linestyle="none",
                marker="o",
                color=color,
                markerfacecolor=color,
                markersize=4.5,
                label=group,
            )
        )

    # Gene labels
    missing_labels = SHOW_LABELS - set(data["gene_symbol"])
    if missing_labels:
        warnings.warn(f"Requested labels absent from plotted data: {sorted(missing_labels)}", stacklevel=2)

    for row in data.iter_rows(named=True):
        if row["gene_symbol"] not in SHOW_LABELS:
            continue

        ax.annotate(
            row["gene_symbol"],
            (row["log2_fc"], row["ctrl_expression"]),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=FONT_SIZE,
            color=ManuscriptColors.BLACK,
        )

    # Significance legend
    handles.extend([
        Line2D(
            [], [],
            linestyle="none",
            marker="o",
            color=ManuscriptColors.BLACK,
            markerfacecolor=ManuscriptColors.BLACK,
            markersize=4.5,
            label=f"Raw p < {P_VAL_THRESHOLD:g}",
        ),
        Line2D(
            [], [],
            linestyle="none",
            marker="o",
            color=ManuscriptColors.BLACK,
            markerfacecolor=ManuscriptColors.WHITE,
            markersize=4.5,
            label=f"Raw p ≥ {P_VAL_THRESHOLD:g} or missing",
        ),
    ])

    # Axes
    ax.set_xlabel(
        r"log$_2$ fold change (Ctrl vs KO)",
        fontsize=X_AXIS_SIZE,
        color=ManuscriptColors.BLACK,
    )
    ax.set_ylabel(
        r"Mean Ctrl expression (log$_2$ intensity)",
        fontsize=FONT_SIZE,
        color=ManuscriptColors.BLACK,
    )

    ax.tick_params(labelsize=SCALE_SIZE, colors=ManuscriptColors.BLACK)
    ax.margins(x=0.15, y=0.12)
    ax.set_facecolor(ManuscriptColors.WHITE)

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    for side in ("bottom", "left"):
        ax.spines[side].set_color(ManuscriptColors.BLACK)

    return handles


def _save_legend(handles: list[Line2D], output_file: Path) -> Path:
    """Save tissue colors and significance markers separately."""

    fig = create_figure(
        width_mm=95,
        height_mm=8 + 5 * len(handles),
        dpi=DPI,
    )

    try:
        fig.legend(
            handles=handles,
            loc="center",
            frameon=False,
            fontsize=FONT_SIZE,
            labelcolor=ManuscriptColors.BLACK,
            handletextpad=0.4,
        )
        return save_svg(fig, output_file)
    finally:
        plt.close(fig)


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def create_aps1_antigens(layout: PanelLayout, output_file: Path) -> Path:
    """Generate the APS-1 antigen panel and separate legend SVG."""

    fig = create_figure(width_mm=layout.width, height_mm=layout.height, dpi=DPI)

    try:
        fig.set_layout_engine("constrained")
        ax = fig.subplots()

        handles = _draw_panel(ax)
        add_panel_label(ax, layout.name)

        legend_file = output_file.with_name(f"{output_file.stem}_legend.svg")
        _save_legend(handles, legend_file)

        return save_svg(fig, output_file)
    finally:
        plt.close(fig)
