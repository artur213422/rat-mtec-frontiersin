from pathlib import Path

import polars as pl
from matplotlib import pyplot as plt

from src.common.feratures.visualization.heatmap.heatmap import create_heatmap
from src.common.feratures.visualization.heatmap.heatmap_models import HeatmapData, HeatmapSpec, HeatmapStyle
from src.common.paths import resolve_dir_from_root
from src.common.types.enums import Colors

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

output_file: Path = resolve_dir_from_root("data/mtec_low/analysis/kallikrein_expression_pattern/results") / "kallikreins.png"
analysed_file: Path = resolve_dir_from_root("data/mtec_low/analysis/kallikrein_expression_pattern/kallikrein_lists") / "kallikreins.csv"

analysed_tsa_df: pl.DataFrame = pl.read_csv(
    analysed_file,
    schema_overrides={
        "ortholog_id": pl.String,
    },
)

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

heatmap_data: HeatmapData = HeatmapData.from_long_frame(
    data=analysed_tsa_df,
    row_label_col="gene_symbol",
    col_label_col="dataset_id",
    value_col="log2_fold_change",
)

heatmap_data.move_col_to_left("rat_mtec_low")
heatmap_data.sort_rows_by_col("rat_mtec_low")

style: HeatmapStyle = HeatmapStyle(
    high_color=Colors.RED,
    low_color=Colors.GREEN,
    center_color=Colors.WHITE,
    nan_cell_color=Colors.BLACK,
    vcenter=0,
    # vmax=3,
    # vmin=-7,
)

spec: HeatmapSpec = HeatmapSpec(
    data=heatmap_data,
    style=style,
)

fig, ax = plt.subplots(
    figsize=(5, 14),
)
fig.suptitle("Kallikreins")

create_heatmap(ax, spec)

##########################################################################################
#                                           SAVE
##########################################################################################

fig.tight_layout()

fig.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight",
)

plt.close(fig)