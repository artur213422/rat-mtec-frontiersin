from itertools import combinations
from pathlib import Path

from src.common.feratures.visualization.venn_diagram.venn_diagram_types import VennDiagramData, VennDiagramSet


def write_overlapped(data: VennDiagramData, file: Path) -> None:
    lines: list[str] = []

    number_of_sets: int = len(data.sets)

    for overlap_size in range(2, number_of_sets + 1):
        for selected_sets in combinations(data.sets, overlap_size):
            overlap: frozenset[str] = _calculate_exact_overlap(selected_sets=selected_sets, all_sets=data.sets,)

            labels: str = " ∩ ".join(
                venn_set.label
                for venn_set in selected_sets
            )

            genes: str = (
                ", ".join(sorted(overlap))
                if overlap
                else "None"
            )

            lines.append(
                f"{labels} ({len(overlap)} genes)"
            )
            lines.append(genes)
            lines.append("")

    file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with file.open(
            mode="w",
            encoding="utf-8",
            newline="\n",
    ) as output:
        output.write("\n".join(lines))


def _calculate_exact_overlap(selected_sets: tuple[VennDiagramSet, ...], all_sets: tuple[VennDiagramSet, ...]) -> frozenset[str]:
    overlap: frozenset[str] = frozenset.intersection(
        *(
            venn_set.data
            for venn_set in selected_sets
        )
    )

    selected_ids: set[int] = {
        id(venn_set)
        for venn_set in selected_sets
    }

    excluded_data: frozenset[str] = frozenset().union(
        *(
            venn_set.data
            for venn_set in all_sets
            if id(venn_set) not in selected_ids
        )
    )

    return overlap - excluded_data
