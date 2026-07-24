from matplotlib.axes import Axes
from matplotlib_venn import venn2, venn3

from src.common.feratures.visualization.venn_diagram.venn_diagram_types import VennDiagramSpec, VennDiagramSet


def create_venn_diagram(ax: Axes,spec: VennDiagramSpec) -> None:
    venn_sets: tuple[VennDiagramSet, ...] = spec.data.sets

    labels: tuple[str, ...] = tuple(
        venn_set.label
        for venn_set in venn_sets
    )

    colors: tuple[str, ...] = tuple(
        venn_set.color
        for venn_set in venn_sets
    )

    values: tuple[set[str], ...] = tuple(
        set(venn_set.data)
        for venn_set in venn_sets
    )

    match len(venn_sets):
        case 2:
            venn2(
                subsets=(values[0], values[1]),
                set_labels=(labels[0], labels[1]),
                set_colors=(colors[0], colors[1]),
                alpha=spec.style.alpha,
                ax=ax,
            )

        case 3:
            venn3(
                subsets=(values[0], values[1], values[2]),
                set_labels=(labels[0], labels[1], labels[2]),
                set_colors=(colors[0], colors[1], colors[2]),
                alpha=spec.style.alpha,
                ax=ax,
            )
        case _:
            raise ValueError(
                "A Venn diagram requires two or three sets."
            )