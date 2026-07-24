from pathlib import Path
import polars as pl
from matplotlib import pyplot as plt

from src.common.feratures.visualization.heatmap.heatmap import create_heatmap
from src.common.feratures.visualization.heatmap.heatmap_models import HeatmapStyle, HeatmapSpec, HeatmapData
from src.common.models.dataset_metadata.model import DatasetCatalog
from src.common.models.dataset_metadata.ops import build_dataset_catalog
from src.common.paths import try_find_file_in_project, resolve_dir_from_root
from src.common.types.enums import Colors
from src.resource_processing.define_tsa_pipeline.internal.tissue_extraction import extract_list_from_json

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


output_folder: Path = resolve_dir_from_root("data/analysis/tsa_expression_pattern/figures")

analysed_tsa_file: Path = try_find_file_in_project("analysed_tsa.csv")
picked_tissues_file: Path = try_find_file_in_project("picked_tissues.json")

analysed_tsa_df: pl.DataFrame = pl.read_csv(analysed_tsa_file)

picked_tissues: list[str] = extract_list_from_json(picked_tissues_file)

datasets: DatasetCatalog = build_dataset_catalog()
rat_dataset: str = datasets.get_rat_dataset().dataset_id

##########################################################################################
#                                       PROCESS DATA
##########################################################################################


for tissue in picked_tissues:


    selected_data: pl.DataFrame = analysed_tsa_df.filter(
        (pl.col("tissue") == tissue) & (pl.col("found_in_dataset") == True)
    )

    print(f"Processing: {tissue}")

    if selected_data.is_empty():
        print(f"{tissue}: 0 rows found. Skipping")
        continue


    heatmap_data: HeatmapData = HeatmapData.from_long_frame(
        data=selected_data,
        row_label_col="gene_symbol",
        col_label_col="dataset",
        value_col="log2_fold_change",
    )

    heatmap_data.move_col_to_left(rat_dataset)
    heatmap_data.sort_rows_by_col(rat_dataset)

    style: HeatmapStyle = HeatmapStyle(
        high_color=Colors.RED,
        low_color=Colors.GREEN,
        center_color=Colors.WHITE,
        nan_cell_color=Colors.BLACK,
        vcenter=0,
        #vmax=3,
        #vmin=-7,

    )

    spec: HeatmapSpec = HeatmapSpec(
        data=heatmap_data,
        style=style
    )

    fig, ax = plt.subplots()
    fig.suptitle(f"{tissue}")

    ax = create_heatmap(ax, spec)

    ##########################################################################################
    #                                           SAVE
    ##########################################################################################
    output_path: Path = output_folder / f"{tissue}.png"

    fig.tight_layout()

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)




for tissue in picked_tissues:
    continue
    selected_data: pl.DataFrame = analysed_tsa_df.filter(
        (pl.col("tissue") == tissue) & (pl.col("found_in_dataset") == True)
    )

    if selected_data.is_empty() and selected_data.height < 5:
        print(f"{tissue}: 0 rows found. Skipping")
        continue

    for dataset in datasets.datasets:
        dataset_id = dataset.dataset_id

        range_summary: pl.DataFrame = (
            analysed_tsa_df
            .filter(pl.col("found_in_dataset"))
            .group_by(
                ["tissue", "dataset"],
                maintain_order=True,
            )
            .agg(
                pl.len().alias("n"),

                pl.col("ctrl_mean").min().alias("ctrl_min"),
                pl.col("ctrl_mean").max().alias("ctrl_max"),

                pl.col("ko_mean").min().alias("ko_min"),
                pl.col("ko_mean").max().alias("ko_max"),

                pl.col("log2_fold_change").min().alias("log2fc_min"),
                pl.col("log2_fold_change").max().alias("log2fc_max"),
                pl.col("log2_fold_change").median().alias("log2fc_median"),
            )
            .with_columns(
                (
                        pl.col("log2fc_max")
                        - pl.col("log2fc_min")
                ).alias("log2fc_range")
            )
            .sort(["tissue", "dataset"])
        )

        range_summary.write_csv(
            output_folder / "range_summary.txt",
            separator="\t",
        )
