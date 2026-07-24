from matplotlib.axes import Axes

from src.common.feratures.visualization.upset_plot.helpers import _get_legend_label
from src.common.feratures.visualization.upset_plot.upset_plot_types import UpSetPlotSpec


def create_upset_plot(ax: Axes,spec: UpSetPlotSpec) -> None:
    style = spec.style
    intersections = spec.data.intersections

    used_labels: set[str] = set()

    for x_position, intersection in enumerate(intersections):
        bottom: int = 0

        for segment in intersection.segments:
            legend_label: str | None = _get_legend_label(
                label=segment.label,
                used_labels=used_labels,
            )

            ax.bar(
                x=x_position,
                height=segment.count,
                bottom=bottom,
                width=style.bar_width,
                color=segment.color.value,
                label=legend_label,
            )

            bottom += segment.count

        if style.show_bar_counts:
            ax.text(
                x=x_position,
                y=intersection.total_count,
                s=str(intersection.total_count),
                horizontalalignment="center",
                verticalalignment="bottom",
            )

    x_positions: list[int] = list(
        range(len(intersections))
    )

    tuple_labels: list[str] = [
        intersection.tuple_label
        for intersection in intersections
    ]

    ax.set_xticks(
        ticks=x_positions,
        labels=tuple_labels,
        rotation=style.x_tick_rotation,
        horizontalalignment="right",
    )

    ax.set_ylabel("Intersection size")

    ax.grid(
        axis="y",
        alpha=0.25,
    )
    ax.set_axisbelow(True)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if used_labels:
        ax.legend(
            frameon=False,
        )
