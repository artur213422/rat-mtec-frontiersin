import numpy as np
import numpy.typing as npt
import polars as pl

from matplotlib.axes import Axes
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

from src.manuscript.shared.config.image_config import (
    FONT_SIZE,
    SCALE_SIZE,
    X_AXIS_SIZE,
    ManuscriptColors,
)
from src.manuscript.shared.heatmap.heatmap_input import ManuscriptHeatmapInput, RowOrder


def draw_manuscript_heatmap(canvas: Axes, heatmap_input: ManuscriptHeatmapInput) -> None:
    if canvas.has_data():
        raise ValueError("Heatmap requires an empty Axes.")

    required: set[str] = {
        heatmap_input.value_col,
        heatmap_input.row_label_col,
        heatmap_input.dataset_id_col,
        heatmap_input.p_col,
        heatmap_input.fdr_col,
    }

    missing: set[str] = required - set(heatmap_input.data.columns)

    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")

    if not heatmap_input.vmin < 0 < heatmap_input.vmax:
        raise ValueError("Color limits must satisfy vmin < 0 < vmax.")

    if not heatmap_input.col_order:
        raise ValueError("Column order cannot be empty.")

    if len(set(heatmap_input.col_order)) != len(heatmap_input.col_order):
        raise ValueError("Column order contains duplicates.")

    # ------------------------------------------------------------------
    # Prepare data
    # ------------------------------------------------------------------

    data: pl.DataFrame = (
        heatmap_input.data
        .select(
            pl.col(heatmap_input.row_label_col).cast(pl.String).alias("_row"),
            pl.col(heatmap_input.dataset_id_col).cast(pl.String).alias("_col"),
            pl.col(heatmap_input.value_col).cast(pl.Float64).alias("_value"),
            pl.col(heatmap_input.p_col).cast(pl.Float64).alias("_p"),
            pl.col(heatmap_input.fdr_col).cast(pl.Float64).alias("_fdr"),
        )
        .filter(pl.col("_col").is_in(heatmap_input.col_order))
    )

    if data.is_empty():
        raise ValueError("No heatmap data available.")

    if data["_row"].null_count():
        raise ValueError("Row labels cannot contain nulls.")

    duplicates = (
        data.group_by("_row", "_col")
        .len()
        .filter(pl.col("len") > 1)
    )

    if not duplicates.is_empty():
        raise ValueError(f"Duplicate heatmap cells:\n{duplicates}")

    missing_datasets = set(heatmap_input.col_order) - set(data["_col"].to_list())

    if missing_datasets:
        raise ValueError(f"Missing datasets: {sorted(missing_datasets)}")

    rows: list[str] = sorted(
        data["_row"].unique().to_list(),
        key=str.casefold,
    )

    if isinstance(heatmap_input.row_order, list):
        requested_rows = heatmap_input.row_order
        if len(requested_rows) != len(set(requested_rows)):
            raise ValueError("Row order contains duplicates.")
        if set(requested_rows) != set(rows):
            raise ValueError("Explicit row order must contain every data row exactly once.")
        rows = list(requested_rows)

    cols: list[str] = list(heatmap_input.col_order)

    row_indices: dict[str, int] = {
        label: index for index, label in enumerate(rows)
    }

    col_indices: dict[str, int] = {
        label: index for index, label in enumerate(cols)
    }

    shape: tuple[int, int] = (len(rows), len(cols))

    matrix: npt.NDArray[np.float64] = np.full(shape, np.nan)
    p_matrix: npt.NDArray[np.float64] = np.full(shape, np.nan)
    fdr_matrix: npt.NDArray[np.float64] = np.full(shape, np.nan)

    for row_label, col_label, value, p_value, fdr in data.iter_rows():
        i: int = row_indices[row_label]
        j: int = col_indices[col_label]

        matrix[i, j] = float(value) if value is not None else np.nan
        p_matrix[i, j] = float(p_value) if p_value is not None else np.nan
        fdr_matrix[i, j] = float(fdr) if fdr is not None else np.nan

    if not np.isfinite(matrix).any():
        raise ValueError("Heatmap has no finite values.")

    # ------------------------------------------------------------------
    # Row ordering
    # ------------------------------------------------------------------

    if heatmap_input.row_order == RowOrder.VALUE_ASCENDING:
        row_means: list[float] = []

        for i in range(len(rows)):
            finite = matrix[i, np.isfinite(matrix[i])]

            mean = float(finite.mean()) if finite.size else float("inf")
            row_means.append(mean)

        sort_values = matrix[:, 0]

        order = sorted(
            range(len(rows)),
            key=lambda i: (
                not np.isfinite(sort_values[i]),
                sort_values[i] if np.isfinite(sort_values[i]) else float("inf"),
                rows[i].casefold(),
            ),
        )

        rows = [rows[i] for i in order]
        matrix = matrix[order, :]
        p_matrix = p_matrix[order, :]
        fdr_matrix = fdr_matrix[order, :]

    # ------------------------------------------------------------------
    # Heatmap colors
    # ------------------------------------------------------------------

    cmap = LinearSegmentedColormap.from_list(
        "manuscript_heatmap",
        [
            heatmap_input.min_color,
            ManuscriptColors.WHITE,
            heatmap_input.max_color,
        ],
        N=256,
    )

    cmap.set_bad(heatmap_input.na_color)

    norm = TwoSlopeNorm(
        vmin=heatmap_input.vmin,
        vcenter=0,
        vmax=heatmap_input.vmax,
    )

    image = canvas.imshow(
        np.ma.masked_invalid(matrix),
        cmap=cmap,
        norm=norm,
        aspect="auto",
        interpolation="nearest",
    )

    # ------------------------------------------------------------------
    # Cell annotations
    # ------------------------------------------------------------------

    for i, j in np.ndindex(matrix.shape):
        value: float = matrix[i, j]

        if not np.isfinite(value):
            continue

        p: float = p_matrix[i, j]
        fdr: float = fdr_matrix[i, j]

        significance = ""

        if np.isfinite(p) and p < heatmap_input.p_threshold:
            significance += "*"

        if np.isfinite(fdr) and fdr < heatmap_input.fdr_threshold:
            significance += "†"

        if heatmap_input.show_values:
            text = f"{value:+.2f}{significance}"
        else:
            text = significance

        if not text:
            continue

        r, g, b, _ = cmap(norm(value))
        luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b

        text_color = (
            ManuscriptColors.WHITE
            if luminance < 0.5
            else ManuscriptColors.BLACK
        )

        canvas.text(
            j,
            i,
            text,
            ha="center",
            va="center",
            fontsize=FONT_SIZE,
            fontweight="bold" if significance else "normal",
            color=text_color,
        )

    # ------------------------------------------------------------------
    # Labels
    # ------------------------------------------------------------------

    display_cols = [
        heatmap_input.col_labels.get(col, col)
        if heatmap_input.col_labels is not None
        else col
        for col in cols
    ]

    canvas.set_xticks(
        range(len(cols)),
        display_cols,
        fontsize=X_AXIS_SIZE,
        fontweight="bold",
        color=ManuscriptColors.BLACK,
    )

    canvas.xaxis.tick_top()

    canvas.set_yticks(
        range(len(rows)),
        rows,
        fontsize=FONT_SIZE,
        color=ManuscriptColors.BLACK,
    )

    for label in canvas.get_yticklabels():
        label.set_fontstyle(
            "italic" if heatmap_input.italic_row_labels else "normal"
        )

    canvas.tick_params(length=0)

    # ------------------------------------------------------------------
    # Cell borders
    # ------------------------------------------------------------------

    canvas.set_xticks(
        np.arange(-0.5, len(cols), 1),
        minor=True,
    )

    canvas.set_yticks(
        np.arange(-0.5, len(rows), 1),
        minor=True,
    )

    canvas.grid(
        which="minor",
        color=ManuscriptColors.WHITE,
        linewidth=1,
    )

    canvas.tick_params(which="minor", length=0)

    for spine in canvas.spines.values():
        spine.set_visible(False)

    # ------------------------------------------------------------------
    # Significance legend
    # ------------------------------------------------------------------

    canvas.text(
        0,
        -0.03,
        f"* p < {heatmap_input.p_threshold:g}; "
        f"† FDR < {heatmap_input.fdr_threshold:g}",
        transform=canvas.transAxes,
        ha="left",
        va="top",
        fontsize=FONT_SIZE,
        color=ManuscriptColors.BLACK,
    )

    # ------------------------------------------------------------------
    # Colorbar
    # ------------------------------------------------------------------

    colorbar_ax = canvas.inset_axes((1.04, 0.55, 0.04, 0.45))

    colorbar = canvas.figure.colorbar(
        image,
        cax=colorbar_ax,
    )

    colorbar.set_label(
        heatmap_input.colorbar_label,
        fontsize=FONT_SIZE,
        color=ManuscriptColors.BLACK,
    )

    colorbar.ax.tick_params(
        labelsize=SCALE_SIZE,
        colors=ManuscriptColors.BLACK,
    )

    # ------------------------------------------------------------------
    # Missing value legend
    # ------------------------------------------------------------------

    na_ax = canvas.inset_axes((1.04, 0.49, 0.04, 0.04))
    na_ax.set_facecolor(heatmap_input.na_color)
    na_ax.set_xticks([])
    na_ax.set_yticks([])

    na_ax.text(
        1.5,
        0.5,
        "NA",
        transform=na_ax.transAxes,
        ha="left",
        va="center",
        fontsize=FONT_SIZE,
        color=ManuscriptColors.BLACK,
    )
