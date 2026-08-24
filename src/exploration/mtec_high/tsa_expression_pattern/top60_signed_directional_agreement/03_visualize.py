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

output_folder: Path = resolve_dir_from_root("data/mtec_high/analysis/tsa_expression_pattern/signed_directional_agreement_against_rat")

directional_file: Path = output_folder / "signed_directional_agreement.csv"


directional_df: pl.DataFrame = pl.read_csv(directional_file)

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

if directional_df.is_empty():
    raise ValueError("Directional agreement table is empty")

heatmap_data: HeatmapData = HeatmapData.from_long_frame(
    data=directional_df,
    row_label_col="tissue_group",
    col_label_col="dataset",
    value_col="agreement",
)

n_rows: int = len(heatmap_data.row_labels)
n_cols: int = len(heatmap_data.col_labels)

fig_width: float = max(5.0, n_cols * 1.2)
fig_height: float = max(5.0, n_rows * 0.45)

style: HeatmapStyle = HeatmapStyle(
    high_color=Colors.RED,
    low_color=Colors.GREEN,
    center_color=Colors.WHITE,
    nan_cell_color=Colors.BLACK,
    vcenter=0.0,
    vmin=-1.0,
    vmax=1.0,
    show_row_labels=True,
)

spec: HeatmapSpec = HeatmapSpec(
    data=heatmap_data,
    style=style,
)

fig, ax = plt.subplots(figsize=(fig_width, fig_height))
fig.suptitle("Signed directional agreement against rat")

create_heatmap(ax, spec)

##########################################################################################
#                                           SAVE
##########################################################################################

output_path: Path = output_folder / "signed_directional_agreement.png"

fig.tight_layout()

fig.savefig(output_path,  dpi=300,  bbox_inches="tight")

plt.close(fig)