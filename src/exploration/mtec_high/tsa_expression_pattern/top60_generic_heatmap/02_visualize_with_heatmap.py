from pathlib import Path
import polars as pl
from matplotlib import pyplot as plt

from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.feratures.visualization.heatmap.heatmap import create_heatmap
from src.common.feratures.visualization.heatmap.heatmap_models import HeatmapStyle, HeatmapSpec, HeatmapData
from src.common.paths import try_find_file_in_project, resolve_dir_from_root
from src.common.types.enums import Colors

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


output_folder: Path = resolve_dir_from_root("data/mtec_high/analysis/tsa_expression_pattern/figures")

tsa_mtec_expression_file: Path = try_find_file_in_project("tsa_mtec_expression.csv")

datasets: list[DatasetMetadata] = get_mtec_high()


tsa_mtec_expression_df: pl.DataFrame = pl.read_csv(tsa_mtec_expression_file)

tissue_groups: list[str] = (
    tsa_mtec_expression_df
    .get_column("tissue_group")
    .unique()
    .sort()
    .to_list()
)

##########################################################################################
#                                       PROCESS DATA
##########################################################################################
for tissue_group in tissue_groups:
    print(f"Processing: {tissue_group}")

    selected_data: pl.DataFrame = tsa_mtec_expression_df.filter(
        (pl.col("tissue_group") == tissue_group)
        & pl.col("found_in_dataset")
    )

    if selected_data.is_empty():
        print(f"{tissue_group}: 0 rows found. Skipping")
        continue

    heatmap_data: HeatmapData = HeatmapData.from_long_frame(
        data=selected_data,
        row_label_col="gene_symbol",
        col_label_col="dataset",
        value_col="log2_fold_change",
    )

    heatmap_data.move_col_to_left("rat_mtec_high")
    heatmap_data.sort_rows_by_col("rat_mtec_high")

    style: HeatmapStyle = HeatmapStyle(
        high_color=Colors.RED,
        low_color=Colors.GREEN,
        center_color=Colors.WHITE,
        nan_cell_color=Colors.BLACK,
        vcenter=0,
        # vmax=3,
        # vmin=-3,
    )

    spec: HeatmapSpec = HeatmapSpec(data=heatmap_data, style=style)

    fig, ax = plt.subplots(figsize=(8, 15))

    fig.suptitle(tissue_group)

    create_heatmap(ax,  spec)

    ##########################################################################################
    #                                           SAVE
    ##########################################################################################
    output_path: Path = output_folder / f"{tissue_group}.png"

    fig.tight_layout()

    fig.savefig(output_path, dpi=300, bbox_inches="tight",)

    plt.close(fig)

