import numpy.typing as npt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.image import AxesImage
import polars as pl

from src.common.feratures.visualization.heatmap.types.heatmap_data import HeatmapData
from src.common.types.enums import Colors
from collections.abc import Mapping
from typing import Literal


class Heatmap:
    heatmap_ax: Axes
    data: HeatmapData

    matrix: np.ndarray[np.float64]
    row_labels: list[str]
    col_labels: list[str]

    colorbar_size: float

    p_threshold: float | None = None
    fdr_threshold: float | None = None
    value_threshold: float | None = None

    colorbar_ticks: npt.NDArray[np.float64] | None

    def __init__(self, heatmap_ax: Axes, data: HeatmapData) -> None:
        if heatmap_ax.has_data():
            raise ValueError("Heatmap requires an empty Axes.")

        self.heatmap_ax: Axes = heatmap_ax
        self.data: HeatmapData = data

        self.colorbar_ticks: npt.NDArray[np.float64] | None = None

        self.matrix = self.data.matrix_for(self.data.value_col)
        self.row_labels: list[str] = self.data.row_labels
        self.col_labels: list[str] = self.data.col_labels

        cmap = LinearSegmentedColormap.from_list(
            "heatmap",
            [
                Colors.DARK_BLUE,
                Colors.OFF_WHITE,
                Colors.DARK_RED,
            ],
            N=256,
        )
        cmap.set_bad(Colors.DARK_GRAY)

        finite_values = self.matrix[np.isfinite(self.matrix)]

        if finite_values.size == 0:
            raise ValueError("Heatmap contains no finite values.")

        vmin = float(finite_values.min())
        vmax = float(finite_values.max())

        self.image: AxesImage = self.heatmap_ax.imshow(
            np.ma.masked_invalid(self.matrix),
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            aspect="auto",
            interpolation="nearest",
        )

        n_rows, n_cols = self.matrix.shape
        self._draw_cell_grid(
            n_rows=n_rows,
            n_cols=n_cols,
        )

        self.heatmap_ax.set_xticks([])
        self.heatmap_ax.set_yticks([])

        for spine in self.heatmap_ax.spines.values():
            spine.set_visible(False)

        self.show_col_labels()
        self.show_row_labels()

    def _draw_cell_grid(self, n_rows: int, n_cols: int,) -> None:
        self.heatmap_ax.set_xticks(
            np.arange(-0.5, n_cols, 1),
            minor=True,
        )

        self.heatmap_ax.set_yticks(
            np.arange(-0.5, n_rows, 1),
            minor=True,
        )

        self.heatmap_ax.grid(
            which="minor",
            color=Colors.WHITE,
            linewidth=2,
        )

        self.heatmap_ax.tick_params(
            which="minor",
            length=0,
        )

    def show_row_labels(self, *,position: Literal["left", "right"] = "left", bold: bool = False,) -> None:
        self.heatmap_ax.set_yticks(
            range(len(self.row_labels)),
            self.row_labels,
            fontweight="bold" if bold else "normal",
        )

        if position == "left":
            self.heatmap_ax.yaxis.tick_left()
        else:
            self.heatmap_ax.yaxis.tick_right()

        self.heatmap_ax.tick_params(axis="y", length=0)

    def show_col_labels(
            self,
            *,
            position: Literal["top", "bottom"] = "top",
            bold: bool = True,
            rotation: float = 0.0,
            labels: Mapping[str, str] | None = None,
    ) -> None:
        display_labels = self.col_labels

        if labels is not None:
            requested_order = list(labels)

            missing = set(self.col_labels) - set(requested_order)
            unknown = set(requested_order) - set(self.col_labels)

            if missing:
                raise ValueError(f"Missing column labels: {sorted(missing)}")

            if unknown:
                raise ValueError(f"Unknown column labels: {sorted(unknown)}")

            indices = [
                self.col_labels.index(label)
                for label in requested_order
            ]

            self.matrix = self.matrix[:, indices]
            self.col_labels = requested_order

            self.image.set_data(
                np.ma.masked_invalid(self.matrix)
            )

            display_labels = [
                labels[label]
                for label in self.col_labels
            ]

        self.heatmap_ax.set_xticks(
            range(len(display_labels)),
            display_labels,
            fontweight="bold" if bold else "normal",
            rotation=rotation,
        )

        if position == "top":
            self.heatmap_ax.xaxis.tick_top()
        else:
            self.heatmap_ax.xaxis.tick_bottom()

        self.heatmap_ax.tick_params(axis="x", length=0)

    def show_significance(
            self,
            value_threshold: float,
            p_value_col: str | None = None,
            fdr_col: str | None = None,
            *,
            p_threshold: float = 0.05,
            fdr_threshold: float = 0.05,
    ) -> None:
        if p_value_col is None and fdr_col is None:
            raise ValueError("At least one significance column must be provided.")

        for significance_col in (p_value_col, fdr_col):
            if significance_col is not None and significance_col not in self.data.data.columns:
                raise ValueError(f"Column {significance_col!r} not found in DataFrame.")

        if value_threshold < 0:
            raise ValueError("value_threshold must be non-negative.")

        self.value_threshold = value_threshold

        self.p_threshold = (
            p_threshold
            if p_value_col is not None
            else None
        )

        self.fdr_threshold = (
            fdr_threshold
            if fdr_col is not None
            else None
        )

        p_matrix: np.ndarray[np.float64] | None = (
            self.data.matrix_for(
                p_value_col,
                col_order=self.col_labels,
            )
            if p_value_col is not None
            else None
        )

        fdr_matrix: np.ndarray[np.float64] | None = (
            self.data.matrix_for(
                fdr_col,
                col_order=self.col_labels,
            )
            if fdr_col is not None
            else None
        )

        both_requested: bool = p_matrix is not None and fdr_matrix is not None


        for row, col in np.ndindex(self.matrix.shape):
            value: float = self.matrix[row, col]

            if not np.isfinite(value) or abs(value) < value_threshold:
                continue

            symbol: str = ""

            if both_requested:
                p_value: float = p_matrix[row, col]
                fdr: float = fdr_matrix[row, col]

                p_significant: bool = np.isfinite(p_value) and p_value < p_threshold
                fdr_significant: bool = np.isfinite(fdr) and fdr < fdr_threshold


                if p_significant and fdr_significant:
                    symbol = r"$\ast\dagger$"

                elif p_significant:
                    symbol = r"$\ast$"

                elif fdr_significant:
                    symbol = r"$\dagger$"


            elif p_matrix is not None:
                p_value = p_matrix[row, col]

                if np.isfinite(p_value) and p_value < p_threshold:
                    symbol = r"$\ast$"

            elif fdr_matrix is not None:
                fdr = fdr_matrix[row, col]

                if np.isfinite(fdr) and fdr < fdr_threshold:
                    symbol = r"$\ast$"

            if symbol:
                self.heatmap_ax.text(
                    col,
                    row,
                    symbol,
                    ha="center",
                    va="center",
                )

    def show_significance_legend(self) -> None:
        if self.p_threshold is None and self.fdr_threshold is None:
            raise RuntimeError("show_significance() or show_value_with_significance() must be called before show_significance_legend().")

        parts: list[str] = []

        if self.p_threshold is not None and self.fdr_threshold is not None:
            parts.append(rf"$\ast$ p < {self.p_threshold:g}")
            parts.append(rf"$\dagger$ FDR < {self.fdr_threshold:g}")

        elif self.p_threshold is not None:
            parts.append(rf"$\ast$ p < {self.p_threshold:g}")

        elif self.fdr_threshold is not None:
            parts.append(rf"$\ast$ FDR < {self.fdr_threshold:g}")

        if self.value_threshold is not None and self.value_threshold > 0:
            parts.append(rf"markers require $|value| \geq {self.value_threshold:g}$")

        self.heatmap_ax.text(
            0.0,
            -0.02,
            "; ".join(parts),
            transform=self.heatmap_ax.transAxes,
            ha="left",
            va="top",
        )


    def show_color_bar(self, title: str, size: float = 0.35, label_distance: int = 5) -> None:
        if not 0 < size <= 1:
            raise ValueError("size must be between 0 and 1.")

        self.colorbar_size = size

        colorbar_ax = self.heatmap_ax.inset_axes(
            (
                1.04,  # x
                1.0 - size,  # y
                0.04,  # width
                size,  # height
            )
        )

        colorbar = self.heatmap_ax.figure.colorbar(self.image, cax=colorbar_ax, ticks=self.colorbar_ticks)

        colorbar.set_label(
            title,
            rotation=90,
            labelpad=label_distance,
            va="bottom",
        )


    def show_na_legend(self) -> None:
        if not hasattr(self, "colorbar_size"):
            raise RuntimeError("show_color_bar() must be called before show_na().")

        na_height = 0.04
        gap = 0.02

        na_ax = self.heatmap_ax.inset_axes(
            (
                1.04,
                1.0 - self.colorbar_size - gap - na_height,
                0.04,
                na_height,
            )
        )

        na_ax.set_facecolor(Colors.DARK_GRAY)

        na_ax.set_xticks([])
        na_ax.set_yticks([])

        for spine in na_ax.spines.values():
            spine.set_visible(True)
            spine.set_color(Colors.BLACK)
            spine.set_linewidth(0.8)

        na_ax.text(
            1.5,
            0.5,
            "NA",
            transform=na_ax.transAxes,
            ha="left",
            va="center",
        )


    def configure_color_scale(self, *, vmin: float, vmax: float, tick: float) -> None:
        if vmin >= vmax:
            raise ValueError("vmin must be smaller than vmax.")

        if tick <= 0:
            raise ValueError("tick must be positive.")

        if not vmin < 0 < vmax:
            raise ValueError("Diverging heatmap scale must contain zero.")

        self.image.set_norm(
            TwoSlopeNorm(
                vmin=vmin,
                vcenter=0.0,
                vmax=vmax,
            )
        )

        self.colorbar_ticks = np.arange(
            vmin,
            vmax + tick / 2,
            tick,
        )

    def show_value(self) -> None:
        for row, col in np.ndindex(self.matrix.shape):
            value: float = self.matrix[row, col]

            if not np.isfinite(value):
                continue

            self.heatmap_ax.text(
                col,
                row,
                f"{value:+.2f}",
                ha="center",
                va="center",
            )

    def show_value_with_significance(
            self,
            value_threshold: float,
            p_value_col: str | None = None,
            fdr_col: str | None = None,
            *,
            p_threshold: float = 0.05,
            fdr_threshold: float = 0.05,
    ) -> None:
        if p_value_col is None and fdr_col is None:
            raise ValueError(
                "At least one significance column must be provided."
            )

        for significance_col in (p_value_col, fdr_col):
            if (
                    significance_col is not None
                    and significance_col not in self.data.data.columns
            ):
                raise ValueError(
                    f"Column {significance_col!r} not found in DataFrame."
                )

        if value_threshold < 0:
            raise ValueError(
                "value_threshold must be non-negative."
            )

        self.value_threshold = value_threshold

        self.p_threshold = (
            p_threshold
            if p_value_col is not None
            else None
        )

        self.fdr_threshold = (
            fdr_threshold
            if fdr_col is not None
            else None
        )

        p_matrix: np.ndarray[np.float64] | None = (
            self.data.matrix_for(
                p_value_col,
                col_order=self.col_labels,
            )
            if p_value_col is not None
            else None
        )

        fdr_matrix: np.ndarray[np.float64] | None = (
            self.data.matrix_for(
                fdr_col,
                col_order=self.col_labels,
            )
            if fdr_col is not None
            else None
        )

        both_requested: bool = (
                p_matrix is not None
                and fdr_matrix is not None
        )

        for row, col in np.ndindex(self.matrix.shape):
            value: float = self.matrix[row, col]

            if not np.isfinite(value):
                continue

            significance: str = ""

            if abs(value) >= value_threshold:

                if both_requested:
                    p_value: float = p_matrix[row, col]
                    fdr: float = fdr_matrix[row, col]

                    if np.isfinite(fdr) and fdr < fdr_threshold:
                        significance = r"$^{\dagger}$"

                    elif np.isfinite(p_value) and p_value < p_threshold:
                        significance = r"$^{\ast}$"

                elif p_matrix is not None:
                    p_value = p_matrix[row, col]

                    if np.isfinite(p_value) and p_value < p_threshold:
                        significance = r"$^{\ast}$"

                elif fdr_matrix is not None:
                    fdr = fdr_matrix[row, col]

                    if np.isfinite(fdr) and fdr < fdr_threshold:
                        significance = r"$^{\ast}$"

            self.heatmap_ax.text(
                col,
                row,
                f"{value:+.2f}{significance}",
                ha="center",
                va="center",
            )