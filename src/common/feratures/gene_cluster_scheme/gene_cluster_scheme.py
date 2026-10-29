from __future__ import annotations


import polars as pl
from matplotlib.axes import Axes
from matplotlib.collections import PatchCollection
from math import isfinite
from matplotlib.colors import to_rgba

from ._data import _GeneClusterDataMixin
from ._renderer import _ClusterRenderingMixin

class GeneClusterScheme(
    _GeneClusterDataMixin,
    _ClusterRenderingMixin
):
    data: pl.DataFrame
    cluster_genes: list[str]
    ax: Axes

    value_col: str
    gene_col: str
    chromosome_col: str
    start_col: str
    stop_col: str

    chromosome: str | None
    chromosomes: list[str]

    gene_collection: PatchCollection

    gap_break_threshold_kb: int
    significance_col: str

    def __init__(
        self,
        ax: Axes,
        data: pl.DataFrame,
        cluster_genes: list[str],
        *,
        value_col: str = "log2_fold_change",
        gene_col: str = "gene_symbol",
        chromosome_col: str = "chromosome",
        start_col: str = "start",
        stop_col: str = "stop",
        gap_break_threshold_kb: int = 500,
        significance_col: str = "p_val"
    ) -> None:
        if ax.has_data() or ax.child_axes:
            raise ValueError("GeneClusterScheme requires an empty Axes.")

        if not cluster_genes:
            raise ValueError("cluster_genes must contain at least one gene.")

        if gap_break_threshold_kb < 0:
            raise ValueError("gap_break_threshold_kb must be nonnegative.")

        self.gap_break_threshold_kb = gap_break_threshold_kb
        self.ax = ax
        self.cluster_genes = cluster_genes

        self.value_col = value_col
        self.gene_col = gene_col
        self.chromosome_col = chromosome_col
        self.start_col = start_col
        self.stop_col = stop_col
        self.significance_col = significance_col

        self._validate_columns(data)

        self.data = self._prepare_data(data)

        self._validate_cluster()

        self.chromosomes = self.data.get_column(self.chromosome_col).unique(maintain_order=True).to_list()
        self.chromosome = self.chromosomes[0] if len(self.chromosomes) == 1 else None

        self._draw()



    def show_title(self, title: str) -> None:
        self.ax.set_title(title)

    def show_color_bar(
        self,
        *,
        vmax: float | None = None,
        vmin: float | None = None,
        vmax_color: str | None = None,
        vmin_color: str | None = None,
        center_color: str | None = None,
    ) -> None:
        self._show_color_bar(vmax, vmin, vmax_color, vmin_color, center_color)

    def show_scatter_plot(self, *, threshold_line: float, highlight_significant: bool = True) -> None:
        """Show zero and both +/- abs(threshold_line) on every chromosome.

        Both thresholds remain visible even when no values cross them.
        A zero threshold leaves only the solid zero reference line.
        """
        self._show_scatter_plot(threshold_line, highlight_significant)

    def configure_gene_border(self, *, thickness: float, color: str) -> None:
        edgecolor = to_rgba(color)
        for collection in self.gene_collections:
            collection.set_linewidth(thickness)
            collection.set_edgecolor(edgecolor)
        self.ax.figure.canvas.draw_idle()
