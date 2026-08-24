from enum import Enum

from src.common.feratures.analysers.deg_extractor import DatasetDegs
from src.common.feratures.visualization.upset_plot.upset_plot_types import UpSetBarSegment

class DegSelection(Enum):
    UPREGULATED = "Upregulated"
    DOWNREGULATED = "Downregulated"
    ALL = "All DEGs"


def get_segments(selection: DegSelection, data: list[DatasetDegs]) -> list[UpSetBarSegment]:
    result: list[UpSetBarSegment] = []

    for dataset in data:
        count: int = 0

        match selection:
            case DegSelection.UPREGULATED:
                count: int = len(dataset.upregulated)

            case DegSelection.DOWNREGULATED:
                count = len(dataset.downregulated)

            case DegSelection.ALL:
                count = len(dataset.upregulated) + len(dataset.downregulated)

        result.append(
            UpSetBarSegment(
                count=count,
                color=dataset.color,
                label=dataset.dataset_id,
            )
        )

    return result