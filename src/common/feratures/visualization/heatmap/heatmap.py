from matplotlib.axes import Axes
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.image import AxesImage

from src.common.feratures.visualization.heatmap.heatmap_models import HeatmapSpec

def create_scalebar(heatmap_ax: Axes, image: AxesImage, label: str) -> None:
    colorbar = heatmap_ax.figure.colorbar(
        image,
        ax=heatmap_ax,
    )

    colorbar.set_label(label)

def create_heatmap(ax: Axes, spec: HeatmapSpec) -> None:
    cmap = LinearSegmentedColormap.from_list(
        "heatmap",
        [
            spec.style.low_color.value,
            spec.style.center_color.value,
            spec.style.high_color.value,
        ],
    )

    if spec.style.nan_cell_color is not None:
        cmap = cmap.with_extremes(
            bad=str(spec.style.nan_cell_color)
        )

    norm = TwoSlopeNorm(
        vmin=spec.style.vmin,
        vcenter=spec.style.vcenter,
        vmax=spec.style.vmax,
    )

    image: AxesImage = ax.imshow(
        spec.data.matrix,
        cmap=cmap,
        norm=norm,
        aspect="auto",
        interpolation="nearest",
    )

    create_scalebar(
        heatmap_ax=ax,
        image=image,
        label="log2 fold change",
    )

    ax.set_xticks(range(len(spec.data.col_labels)))
    ax.set_xticklabels(
        spec.data.col_labels,
        rotation=90,
    )

    if spec.style.show_row_labels:
        ax.set_yticks(range(len(spec.data.row_labels)))
        ax.set_yticklabels(spec.data.row_labels)
