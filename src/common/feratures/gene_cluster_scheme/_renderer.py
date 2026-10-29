from __future__ import annotations

import matplotlib as mpl
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter, MaxNLocator
from matplotlib.transforms import Bbox, TransformedBbox


class _ClusterRenderingMixin:
    """Rendering and shared color-scale implementation for GeneClusterScheme."""

    def _draw(self) -> None:
        values = self.data.get_column(self.value_col).to_numpy().astype(float)
        finite = values[np.isfinite(values)]
        if not finite.size:
            raise ValueError(f"{self.value_col!r} contains no finite values.")
        absolute_max = float(np.max(np.abs(finite))) or 1.0
        self.cmap = mpl.colormaps["RdBu_r"].copy()
        self.cmap.set_bad("lightgray")
        self.norm = TwoSlopeNorm(vmin=-absolute_max, vcenter=0, vmax=absolute_max)
        self.color_bar = None
        rows = self._chromosome_segments()
        # Equal bp-to-display scale for all segments, independent of their origin.
        spans = [[float(part["_plot_stop"].max() - part["_plot_start"].min())
                  for part in parts] for _, parts in rows]
        if any(sum(widths) <= 0 for widths in spans):
            raise ValueError("Cluster genomic width must be positive.")
        padding = max(max(map(sum, spans)) * 0.03, 10_000.0)
        gap_width = padding * 2
        total_width = max(sum(widths) + 2 * padding * len(widths)
                          + gap_width * (len(widths) - 1) for widths in spans)
        self.ax.set_axis_off()
        self.chromosome_axes = []
        self.segment_axes = []
        self.segment_bounds = []
        self.gene_collections = []
        self.scatter_axes = []
        self.scatter_collections = []
        height = 1 / len(rows)
        for index, ((chromosome, parts), widths) in enumerate(zip(rows, spans)):
            row_ax = self.ax.inset_axes([0, 1 - (index + 1) * height, 1, height])
            row_ax.set_axis_off()
            row_ax.text(0, 0.98, str(chromosome), transform=row_ax.transAxes,
                        ha="left", va="top", fontsize=10)
            self.chromosome_axes.append(row_ax)
            axes = []
            row_bounds = []
            left = 0.0
            for part, width in zip(parts, widths):
                fraction = (width + 2 * padding) / total_width
                segment_ax = row_ax.inset_axes([left, 0.27, fraction, 0.53])
                row_bounds.append((left, fraction))
                self._draw_segment(segment_ax, part, padding)
                axes.append(segment_ax)
                left += fraction
                if len(axes) < len(parts):
                    # Dashed connector explicitly denotes omitted genomic distance.
                    baseline = 0.27 + 0.53 * (0.15 / 1.5)
                    row_ax.plot([left - padding / total_width,
                                 left + (gap_width + padding) / total_width],
                                [baseline, baseline], transform=row_ax.transAxes,
                                color="black", linestyle="--", linewidth=1.2)
                    left += gap_width / total_width
            row_ax.text(0.5, 0.01, "Genomic position (Mb)",
                        transform=row_ax.transAxes, ha="center", va="bottom", fontsize=9)
            self.segment_axes.append(axes)
            self.segment_bounds.append(row_bounds)
        # Preserve the existing public mappable for callers using fig.colorbar.
        self.gene_collection = self.gene_collections[0]

    def _draw_segment(self, ax, data, padding) -> None:
        start = float(data["_plot_start"].min())
        stop = float(data["_plot_stop"].max())
        ax.hlines(0, start, stop, color="black", linewidth=2, zorder=1)
        rectangles = [Rectangle((float(row["_plot_start"]), 0.15),
                                float(row["_plot_stop"] - row["_plot_start"]), 0.5)
                      for row in data.iter_rows(named=True)]
        collection = PatchCollection(rectangles, cmap=self.cmap, norm=self.norm,
                                     edgecolor="black", linewidth=0.8, zorder=3)
        collection.set_array(np.ma.masked_invalid(
            data.get_column(self.value_col).to_numpy().astype(float)))
        ax.add_collection(collection)
        self.gene_collections.append(collection)
        for row in data.iter_rows(named=True):
            ax.text(float(row["_midpoint"]), 0.73, str(row[self.gene_col]),
                    ha="left", va="center", rotation=90, rotation_mode="anchor",
                    fontsize=8, zorder=4)
        ax.set_xlim(start - padding, stop + padding)
        ax.set_ylim(-0.15, 1.35)
        ax.set_yticks([])
        ax.xaxis.set_major_locator(MaxNLocator(nbins=3, min_n_ticks=2))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x / 1_000_000:g}"))
        ax.tick_params(axis="x", labelsize=8)
        for side in ("left", "right", "top"):
            ax.spines[side].set_visible(False)

    def _show_scatter_plot(self, threshold_line, highlight_significant) -> None:
        if not np.isfinite(threshold_line):
            raise ValueError("threshold_line must be finite.")
        if highlight_significant and self.significance_col not in self.data.columns:
            raise ValueError(f"Missing significance column: {self.significance_col!r}")
        # Validate before replacing existing artists.
        if highlight_significant:
            self.data[self.significance_col].to_numpy().astype(float)
        values = self.data[self.value_col].to_numpy().astype(float)
        finite = values[np.isfinite(values)]
        threshold = abs(float(threshold_line))
        threshold_width = 0.8
        # Always include zero and both thresholds, even without crossings.
        low = min(float(finite.min()), -threshold)
        high = max(float(finite.max()), threshold)
        thresholds = [-threshold, threshold] if threshold > 0 else []
        margin = (high - low) * 0.1 or max(abs(low) * 0.1, 0.1)
        for ax in self.scatter_axes:
            ax.remove()
        self.scatter_axes = []
        self.scatter_collections = []
        for row_ax, segment_axes, row_bounds, (_, parts) in zip(
                self.chromosome_axes, self.segment_axes, self.segment_bounds,
                self._chromosome_segments()):
            for segment_ax, (left, width), part in zip(segment_axes, row_bounds, parts):
                # Use the same horizontal bounds and shared x-axis, including breaks.
                scheme_bounds = Bbox.from_bounds(left, 0.27, width, 0.35)
                segment_ax.set_axes_locator(
                    lambda ax, renderer, bounds=scheme_bounds, parent=row_ax:
                    TransformedBbox(bounds, parent.transAxes +
                                    ax.figure.transSubfigure.inverted()))
                # Increase the gap above the scheme from 0.10 to 0.15 row units.
                scatter_ax = row_ax.inset_axes([left, 0.77, width, 0.22], sharex=segment_ax)
                self.scatter_axes.append(scatter_ax)
                # Overlapping genes can have midpoints in a different order
                # from their starts. Connect dots in plotted x order, keeping
                # values and significance attached to the same gene.
                part = part.sort("_midpoint", maintain_order=True)
                x = part["_midpoint"].to_numpy().astype(float)
                y = part[self.value_col].to_numpy().astype(float)
                y[~np.isfinite(y)] = np.nan
                scatter_ax.plot(x, y, color="gray", linewidth=0.8, zorder=1)
                significant = np.ones(len(y), dtype=bool)
                if highlight_significant:
                    p = part[self.significance_col].to_numpy().astype(float)
                    significant = np.isfinite(p) & (p < 0.05)
                points = scatter_ax.scatter(x[significant], y[significant],
                                            c=y[significant], cmap=self.cmap,
                                            norm=self.norm, s=16, zorder=3)
                self.scatter_collections.append(points)
                scatter_ax.scatter(x[~significant], y[~significant],
                                   color="gray", s=16, zorder=3)
                scatter_ax.axhline(0, color="#555555", linestyle="-",
                                   linewidth=threshold_width * 0.5, zorder=2)
                for level in thresholds:
                    scatter_ax.axhline(level, color="red", linestyle="--",
                                       linewidth=threshold_width, zorder=2)
                scatter_ax.set_ylim(low - margin, high + margin)
                scatter_ax.tick_params(axis="both", labelsize=8)
                scatter_ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
                # The shared x-axis inherits the scheme's Mb tick formatter.
                scatter_ax.set_xlabel("Genomic position (Mb)", fontsize=8)
                if segment_ax is segment_axes[0]:
                    scatter_ax.set_ylabel(self.value_col, fontsize=8)
                for side in ("right", "top"):
                    scatter_ax.spines[side].set_visible(False)
        self.ax.figure.canvas.draw_idle()

    def _show_color_bar(self, vmax, vmin, vmax_color, vmin_color, center_color) -> None:
        norm = TwoSlopeNorm(vmin=self.norm.vmin if vmin is None else vmin,
                            vcenter=0,
                            vmax=self.norm.vmax if vmax is None else vmax)
        cmap = self.cmap
        if any(color is not None for color in (vmin_color, center_color, vmax_color)):
            cmap = LinearSegmentedColormap.from_list("gene_cluster", [
                vmin_color if vmin_color is not None else cmap(0.0),
                center_color if center_color is not None else cmap(0.5),
                vmax_color if vmax_color is not None else cmap(1.0)])
            cmap.set_bad("lightgray")
        self.norm, self.cmap = norm, cmap
        for collection in self.gene_collections + self.scatter_collections:
            collection.set_norm(norm)
            collection.set_cmap(cmap)
        if self.color_bar is None:
            self.color_bar = self.ax.figure.colorbar(
                self.gene_collection, ax=self.ax, label="Log2 fold change", pad=0.02)
        else:
            self.color_bar.update_normal(self.gene_collection)
