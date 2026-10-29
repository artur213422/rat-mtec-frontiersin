"""Rat genomic-cluster expression plots ordered by gene start coordinate."""
from pathlib import Path

import numpy as np
import polars as pl
from matplotlib.axes import Axes
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.figure import Figure

from src.common.dataset_catalog.api.get_datasets import get_rat_tec_high, read_file_with_ortholog_id
from src.common.paths import PROJECT_ROOT
from src.manuscript.figures.layout.panel_layout import PanelLayout
from src.manuscript.shared.config.image_config import FONT_SIZE, SCALE_SIZE, TITLE_SIZE, X_AXIS_SIZE, ManuscriptColors, DPI
from src.manuscript.shared.fig.fig import create_figure, save_svg

clusters_file = PROJECT_ROOT / "data/mtec_high/derived_data/genomics/window_sliding/standardized/clusters.csv"

show_clusters = {
    "Mast cell proteases": "Mcpt1l1|Mcpt10|Mcpt8|Mcpt4|Mcpt4l1|Mcpt1|RGD1562035|Mcpt1l3|Mcpt9|Mcpt2",
    "Kallikreins": "Klk5|Klk4|LOC687880|LOC102546874|Klk1c6|Klk1c2|Klk1c9|Klk1|LOC100363199|Klk1c3",
    "Antiviral OAS": "Oas1a|Oas1f|Oas1b|Oas1e|Oas1k|Oas3|Oas2",
    "Insulin locus": "Igf2|LOC102548568|Ins2|Th",
    "Defensins": "Defa5|Defa24|Defal1|LOC102555283|LOC100359767|Defa6|Np4|LOC100359819|Defa8|Defa9",
    "GBP (ISGs)": "Gbp5|LOC685067|Gbp4|Gbp2|Gbp1|Gbp3|LOC100911495",
    "Olfactory receptors": "Olr1553|Olr1555|Olr1557|Olr1558|Olr1559|Olr1560|Olr1561|Olr1562|Cldnd1|Gpr15",
    "Steroid reductases": "Akr1c1|Akr1c3l1|Akr1c12l1|Akr1c12|Akr1c19|Akr1c13|Akr1c3|Akr1c14|Ucn3|Tubal3",
}


DATASET_ID: str = "rat_mtec_high"

SIGNIFICANCE_THRESHOLD = 0.05
FOLD_CHANGE_THRESHOLD = 1.0
GAP = 0.045
MARKER_SIZE = 20


def _plot_cluster(
    ax: Axes, cluster_name: str, genes: list[str], cluster_df: pl.DataFrame,
    *, dataset_id: str = DATASET_ID,
) -> None:
    """Draw one genomic cluster expression scatter plot."""

    data = cluster_df.filter(
        (pl.col("dataset_id") == dataset_id)
        & pl.col("gene_symbol").is_in(genes)
    )

    if data.is_empty():
        raise ValueError(f"No genes found for cluster {cluster_name}")

    missing = set(genes) - set(data["gene_symbol"].to_list())
    if missing:
        raise ValueError(f"Missing genes in cluster {cluster_name}: {sorted(missing)}")
    if data.height != len(genes):
        raise ValueError(f"Duplicate gene rows in cluster {cluster_name}")
    if data["chromosome"].n_unique() != 1:
        raise ValueError(f"Genes in {cluster_name} do not belong to one chromosome")
    if "window_id" in data.columns and data["window_id"].n_unique() != 1:
        raise ValueError(f"Genes in {cluster_name} do not belong to one genomic window")
    if data.filter(
        pl.col("start").is_null()
        | pl.col("stop").is_null()
        | (pl.col("start") > pl.col("stop"))
    ).height:
        raise ValueError(f"Missing or invalid gene coordinates in {cluster_name}")

    # Window boundaries describe the entire cluster, not individual genes.
    # Use starts rather than midpoints so overlapping gene spans retain their
    # genomic start order (for example, the broad Mcpt1l1 annotation).
    data = data.sort(["start", "stop", "gene_symbol"]).with_columns(
        (pl.col("start") / 1_000_000).alias("position_mb")
    )

    x = data["position_mb"].to_numpy()
    y = data["log2_fold_change"].to_numpy()
    p = data["p_val"].to_numpy()

    significant = (
        np.isfinite(p)
        & (p < SIGNIFICANCE_THRESHOLD)
    )

    limit = max(
        float(np.abs(y).max()),
        FOLD_CHANGE_THRESHOLD,
    )

    norm = TwoSlopeNorm(
        vmin=-limit,
        vcenter=0,
        vmax=limit,
    )

    cmap = LinearSegmentedColormap.from_list(
        "cluster_fc",
        [
            ManuscriptColors.BLUE,
            ManuscriptColors.WHITE,
            ManuscriptColors.RED,
        ],
    )

    # connecting line follows genomic order
    ax.plot(
        x,
        y,
        color=ManuscriptColors.BLACK,
        alpha=0.3,
        linewidth=0.8,
        zorder=1,
    )

    ax.scatter(
        x[significant],
        y[significant],
        c=y[significant],
        cmap=cmap,
        norm=norm,
        s=MARKER_SIZE,
        zorder=3,
    )

    ax.scatter(
        x[~significant],
        y[~significant],
        color=ManuscriptColors.BLACK,
        alpha=0.4,
        s=MARKER_SIZE,
        zorder=3,
    )

    ax.axhline(
        0,
        color=ManuscriptColors.BLACK,
        alpha=0.4,
        linewidth=0.5,
    )

    ax.axhline(
        FOLD_CHANGE_THRESHOLD,
        color=ManuscriptColors.RED,
        linestyle="--",
        linewidth=0.6,
        alpha=0.6,
    )

    ax.axhline(
        -FOLD_CHANGE_THRESHOLD,
        color=ManuscriptColors.RED,
        linestyle="--",
        linewidth=0.6,
        alpha=0.6,
    )

    for row in data.iter_rows(named=True):
        ax.annotate(
            row["gene_symbol"],
            (
                row["position_mb"],
                row["log2_fold_change"],
            ),
            xytext=(2, 5),
            textcoords="offset points",
            rotation=90,
            ha="left",
            va="bottom",
            fontsize=FONT_SIZE,
            color=ManuscriptColors.BLACK,
        )

    padding = max(
        float(np.ptp(x)) * 0.08,
        0.005,
    )

    ax.set_xlim(
        float(x.min()) - padding,
        float(x.max()) + padding,
    )

    # Positive clusters can include small negative changes near zero.
    # Remove unused space below -1 only when every point remains visible.
    lower_limit = -1.0 if float(y.min()) >= -1.0 and float(y.max()) > 0 else -limit * 1.3
    ax.set_ylim(
        lower_limit,
        limit * 1.3,
    )

    chromosome = data["chromosome"][0]

    ax.set_title(
        f"{cluster_name} · {chromosome}",
        fontsize=TITLE_SIZE,
        color=ManuscriptColors.BLACK,
    )

    ax.set_xlabel(
        "Gene start (Mb)",
        fontsize=FONT_SIZE,
        color=ManuscriptColors.BLACK,
    )

    # Keep the source start coordinates visible without an additive offset.
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)

    ax.set_ylabel(
        "log₂ fold change (KO/WT)",
        fontsize=X_AXIS_SIZE,
        color=ManuscriptColors.BLACK,
    )

    ax.tick_params(
        labelsize=SCALE_SIZE,
        colors=ManuscriptColors.BLACK,
    )

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    ax.spines["left"].set_color(ManuscriptColors.BLACK)

    ax.spines["bottom"].set_color(ManuscriptColors.BLACK)



def create_genomic_clusters(layout: PanelLayout, output_file: Path) -> Path:
    fig: Figure = create_figure(width_mm=layout.width, height_mm=layout.height, dpi=DPI)

    fig.set_layout_engine("constrained")
    axes = fig.subplots(4, 2).flat

    selected_genes = [gene for genes in show_clusters.values() for gene in genes.split("|")]
    coordinates = (
        get_rat_tec_high().get_standardized_df()
        .filter(pl.col("gene_symbol").is_in(selected_genes))
        .select("gene_symbol", "chromosome", "start", "stop")
    )
    # Gene symbols, not ortholog IDs: several rat genes share an ortholog ID
    # and some genes in these clusters have no ortholog ID at all.
    cluster_df = (
        read_file_with_ortholog_id(clusters_file)
        .filter((pl.col("dataset_id") == DATASET_ID) & pl.col("gene_symbol").is_in(selected_genes))
        .join(coordinates, on=["gene_symbol", "chromosome"], how="left", validate="m:1")
    )

    for ax, (cluster_name, genes) in zip(axes, show_clusters.items(),):
        gene_list = genes.split("|")

        _plot_cluster(ax, cluster_name, gene_list, cluster_df)



    save_svg(fig, output_file)

    return output_file
