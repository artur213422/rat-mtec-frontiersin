from matplotlib.axes import Axes
from matplotlib_venn import venn2, venn3

from src.common.feratures.visualization.venn_diagram.venn_diagram_types import VennDiagramSet, VennDiagramSpec



def _format_region_label(region_size: int, union_size: int) -> str:
    percentage = (
        100.0 * region_size / union_size
        if union_size > 0
        else 0.0
    )

    return f"{region_size:,}\n({percentage:.1f}%)"


def create_venn_diagram(ax: Axes, spec: VennDiagramSpec) -> None:
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

    union_size = len(set().union(*values))

    region_values: dict[str, set[str]] = {}

    match len(venn_sets):
        case 2:
            first, second = values

            region_values = {
                "10": first - second,
                "01": second - first,
                "11": first & second,
            }

            venn = venn2(
                subsets=(first, second),
                set_labels=(
                    f"{labels[0]}\n(n={len(first):,})",
                    f"{labels[1]}\n(n={len(second):,})",
                ),
                set_colors=(colors[0], colors[1]),
                alpha=spec.style.alpha,
                ax=ax,
            )

        case 3:
            first, second, third = values

            region_values = {
                "100": first - second - third,
                "010": second - first - third,
                "001": third - first - second,
                "110": (first & second) - third,
                "101": (first & third) - second,
                "011": (second & third) - first,
                "111": first & second & third,
            }

            venn = venn3(
                subsets=(first, second, third),
                set_labels=(
                    f"{labels[0]}\n(n={len(first):,})",
                    f"{labels[1]}\n(n={len(second):,})",
                    f"{labels[2]}\n(n={len(third):,})",
                ),
                set_colors=(colors[0], colors[1], colors[2]),
                alpha=spec.style.alpha,
                ax=ax,
            )

        case _:
            raise ValueError(
                "A Venn diagram requires two or three sets."
            )

    for region_id, region_genes in region_values.items():
        region_label = venn.get_label_by_id(region_id)

        # Some geometrically empty regions do not have a label object.
        if region_label is None:
            continue

        region_label.set_text(
            _format_region_label(
                region_size=len(region_genes),
                union_size=union_size,
            )
        )