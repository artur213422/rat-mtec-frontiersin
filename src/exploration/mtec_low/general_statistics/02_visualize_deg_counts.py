from pathlib import Path

import polars as pl
from matplotlib import pyplot as plt

from src.common.dataset_catalog.api.get_datasets import get_mtec_low
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.feratures.analysers.deg_extractor import get_degs, DatasetDegs, print_deg_summary
from src.common.feratures.visualization.upset_plot.upset_plot import create_upset_plot
from src.common.feratures.visualization.upset_plot.upset_plot_types import UpSetPlotSpec, UpSetPlotStyle, UpSetPlotData, \
    UpSetIntersection, UpSetBarSegment

from src.common.paths import resolve_dir_from_root
from src.common.types.enums import Colors
from src.exploration.mtec_high.general_statistics.utils.colors import get_dataset_color
from src.exploration.mtec_high.general_statistics.utils.get_upset_segments import get_segments, DegSelection

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


output_folder: Path = resolve_dir_from_root("data/mtec_low/analysis/general_statistics")
datasets: list[DatasetMetadata] = get_mtec_low()

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

degs: list[DatasetDegs] = []


for dataset in datasets:
    print("Running dataset: ", dataset)
    tec_df: pl.DataFrame = dataset.get_standardized_df()

    color: Colors = get_dataset_color(dataset.dataset_id)
    deg: DatasetDegs =  get_degs(dataset.dataset_id, tec_df, color)
    degs.append(deg)

print("Running here")
intersections: list[UpSetIntersection] = []

decreased_deg: list[UpSetBarSegment] = get_segments(DegSelection.DOWNREGULATED, degs)
increased_deg: list[UpSetBarSegment] = get_segments(DegSelection.UPREGULATED, degs)
all_deg: list[UpSetBarSegment] = get_segments(DegSelection.ALL, degs)

intersections.append(
    UpSetIntersection(
        tuple_label="Downregulated DEGs",
        segments=tuple(decreased_deg)
    )
)

intersections.append(
    UpSetIntersection(
        tuple_label="Upregulated DEGs",
        segments=tuple(increased_deg)
    )
)

intersections.append(
    UpSetIntersection(
        tuple_label="All DEGs",
        segments=tuple(all_deg)
    )
)


data: UpSetPlotData = UpSetPlotData(
    intersections=tuple(intersections)
)

##########################################################################################
#                                       CREATE IMAGE
##########################################################################################

spec: UpSetPlotSpec = UpSetPlotSpec(
    data=data,
    style=UpSetPlotStyle()
)

fig, ax = plt.subplots(
    figsize=(10.0, 14.0),
)

create_upset_plot(ax, spec)

fig.suptitle("Deg overlap")


##########################################################################################
#                                       SAVE
##########################################################################################

print_deg_summary(degs)

fig.tight_layout()

fig.savefig(
    output_folder / "deg_count.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(fig)