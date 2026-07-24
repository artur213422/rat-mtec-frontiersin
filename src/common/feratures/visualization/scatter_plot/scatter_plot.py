from __future__ import annotations

import numpy as np
import numpy.typing as npt
from matplotlib.axes import Axes

from src.common.feratures.visualization.scatter_plot.scatter_plot_types import (
    CategoricalScatterPlotSpec,
)


def create_categorical_scatter_plot(ax: Axes, spec: CategoricalScatterPlotSpec,) -> None:
    x_positions: npt.NDArray[np.int64] = np.arange(
        len(spec.data.categories),
        dtype=np.int64,
    )

    ax.plot(
        x_positions,
        spec.data.y_values,
        color=spec.style.marker_color.value,
        alpha=spec.style.marker_alpha,
        linewidth=1.5,
    )

    ax.scatter(
        x=x_positions,
        y=spec.data.y_values,
        color=spec.style.marker_color.value,
        s=spec.style.marker_size,
        alpha=spec.style.marker_alpha,
    )

    ax.set_xticks(
        ticks=x_positions,
        labels=spec.data.categories,
    )