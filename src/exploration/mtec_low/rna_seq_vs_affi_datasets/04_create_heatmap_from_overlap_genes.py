from pathlib import Path
import polars as pl
from matplotlib import pyplot as plt

from src.common.feratures.visualization.heatmap.heatmap import create_heatmap
from src.common.feratures.visualization.heatmap.heatmap_models import HeatmapStyle, HeatmapSpec, HeatmapData

from src.common.paths import try_find_file_in_project, resolve_dir_from_root
from src.common.types.enums import Colors

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


output_file: Path = resolve_dir_from_root("data/mtec_low/analysis/affy_vs_rna_seq/heatmap") / "overlapped_genes.png"

analysed_file: Path = try_find_file_in_project("union_genes.csv")

overlapped_genes_df: pl.DataFrame = pl.read_csv(analysed_file)

##########################################################################################
#                                       PROCESS DATA
##########################################################################################




heatmap_data: HeatmapData = HeatmapData.from_long_frame(
    data=overlapped_genes_df,
    row_label_col="gene_symbol",
    col_label_col="dataset_id",
    value_col="log2_fold_change",
)

heatmap_data.sort_rows_by_col("mouse_gse205410_mtec_low")


style: HeatmapStyle = HeatmapStyle(
        high_color=Colors.RED,
        low_color=Colors.GREEN,
        center_color=Colors.WHITE,
        nan_cell_color=Colors.BLACK,
        vcenter=0,
        show_row_labels=False,
        # vmax=3,
        # vmin=-7,
)

spec: HeatmapSpec = HeatmapSpec(
        data=heatmap_data,
        style=style
)

fig, ax = plt.subplots(
    figsize=(4, 50),
)

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


