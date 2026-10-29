from matplotlib.axes import Axes

from src.common.feratures.visualization.upset_plot.helpers import _get_legend_label
from src.common.feratures.visualization.upset_plot.upset_plot_types import UpSetBarSegment, UpSetIntersection, UpSetPlotSpec
from src.common.types.enums import Colors


def _collapse_intersections_by_tuple_label(intersections: tuple[UpSetIntersection, ...]) -> tuple[UpSetIntersection, ...]:
    """
    Combine intersections having the same tuple_label.

    Segments with the same color and legend label are summed.
    Original tuple-label and segment ordering is preserved.
    """

    grouped_segment_counts: dict[
        str,
        dict[
            tuple[Colors, str | None],
            int,
        ],
    ] = {}

    for intersection in intersections:
        if intersection.tuple_label not in grouped_segment_counts:
            grouped_segment_counts[
                intersection.tuple_label
            ] = {}

        segment_counts = grouped_segment_counts[
            intersection.tuple_label
        ]

        for segment in intersection.segments:
            segment_key: tuple[
                Colors,
                str | None,
            ] = (
                segment.color,
                segment.label,
            )

            segment_counts[segment_key] = (
                segment_counts.get(
                    segment_key,
                    0,
                )
                + segment.count
            )

    collapsed_intersections: list[
        UpSetIntersection
    ] = []

    for tuple_label, segment_counts in grouped_segment_counts.items():
        collapsed_segments: tuple[
            UpSetBarSegment,
            ...,
        ] = tuple(
            UpSetBarSegment(
                count=count,
                color=color,
                label=label,
            )
            for (
                color,
                label,
            ), count in segment_counts.items()
        )

        collapsed_intersections.append(
            UpSetIntersection(
                tuple_label=tuple_label,
                segments=collapsed_segments,
            )
        )

    return tuple(collapsed_intersections)


def create_upset_plot(ax: Axes, spec: UpSetPlotSpec) -> None:
    style = spec.style

    intersections: tuple[
        UpSetIntersection,
        ...,
    ] = spec.data.intersections

    if style.collapse_duplicate_tuple_labels:
        intersections = (
            _collapse_intersections_by_tuple_label(
                intersections
            )
        )

    used_labels: set[str] = set()

    for x_position, intersection in enumerate(intersections):
        bottom: int = 0

        for segment in intersection.segments:
            legend_label: str | None = (
                _get_legend_label(
                    label=segment.label,
                    used_labels=used_labels,
                )
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