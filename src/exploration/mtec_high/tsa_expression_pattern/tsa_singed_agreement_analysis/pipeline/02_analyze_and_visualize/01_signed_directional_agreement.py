from pathlib import Path

import polars as pl
from matplotlib import pyplot as plt

from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.feratures.visualization.heatmap.heatmap import create_heatmap
from src.common.feratures.visualization.heatmap.heatmap_models import (
    HeatmapData,
    HeatmapSpec,
    HeatmapStyle,
)
from src.common.paths import resolve_dir_from_root
from src.common.types.enums import Colors
from src.exploration.mtec_high.tsa_expression_pattern.tsa_singed_agreement_analysis.types.api_for_standardized_data import (
    AireStatus,
    get_tsa_changed_in_any_dataset_df,
)

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

tsa_df: pl.DataFrame = get_tsa_changed_in_any_dataset_df()

datasets: list[DatasetMetadata] = get_mtec_high()

output_folder: Path = resolve_dir_from_root("data/mtec_high/analysis/tsa_expression_pattern/tsa_singed_agreement_analysis/results")

##########################################################################################
#                                       SETTINGS
##########################################################################################

rat_dataset: str = "rat_mtec_high"

mouse_datasets: list[str] = [
    dataset.dataset_id
    for dataset in datasets
    if dataset.dataset_id != rat_dataset
]

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

if tsa_df.is_empty():
    raise ValueError("Changed TSA table is empty")

rat_df: pl.DataFrame = (
    tsa_df
    .filter(pl.col("dataset") == rat_dataset)
    .select(
        "tissue_group",
        "gene_symbol",
        pl.col("aire_status_dataset").alias("rat_aire_status"),
    )
)

result_frames: list[pl.DataFrame] = []

for mouse_dataset in mouse_datasets:
    mouse_df: pl.DataFrame = (
        tsa_df
        .filter(pl.col("dataset") == mouse_dataset)
        .select(
            "tissue_group",
            "gene_symbol",
            pl.col("aire_status_dataset").alias("mouse_aire_status"),
        )
    )

    comparison_df: pl.DataFrame = (
        rat_df
        .join(
            mouse_df,
            on=["tissue_group", "gene_symbol"],
            how="inner",
        )
        .with_columns(
            pl.when(
                (
                    (pl.col("rat_aire_status") == AireStatus.DOWNREGULATED)
                    & (pl.col("mouse_aire_status") == AireStatus.DOWNREGULATED)
                )
                |
                (
                    (pl.col("rat_aire_status") == AireStatus.UPREGULATED)
                    & (pl.col("mouse_aire_status") == AireStatus.UPREGULATED)
                )
            )
            .then(1.0)
            .when(
                (
                    (pl.col("rat_aire_status") == AireStatus.DOWNREGULATED)
                    & (pl.col("mouse_aire_status") == AireStatus.UPREGULATED)
                )
                |
                (
                    (pl.col("rat_aire_status") == AireStatus.UPREGULATED)
                    & (pl.col("mouse_aire_status") == AireStatus.DOWNREGULATED)
                )
            )
            .then(-1.0)
            .otherwise(0.0)
            .alias("_directional_agreement")
        )
        .group_by("tissue_group")
        .agg(
            pl.col("_directional_agreement")
            .mean()
            .alias("agreement")
        )
        .with_columns(
            pl.lit(mouse_dataset).alias("dataset")
        )
        .select(
            "tissue_group",
            "dataset",
            "agreement",
        )
    )

    result_frames.append(comparison_df)

if not result_frames:
    raise ValueError("No mouse datasets available for comparison")

result_df: pl.DataFrame = (
    pl.concat(result_frames)
    .sort(["tissue_group", "dataset"])
)

##########################################################################################
#                                       QUICK VIEW
##########################################################################################

print(result_df)

heatmap_data: HeatmapData = HeatmapData.from_long_frame(
    data=result_df,
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

fig, ax = plt.subplots(
    figsize=(fig_width, fig_height),
)

fig.suptitle("Signed directional agreement against rat")

create_heatmap(
    ax,
    spec,
)

fig.tight_layout()

##########################################################################################
#                                           SAVE
##########################################################################################

output_path: Path = output_folder / "signed_directional_agreement_quick_view.png"

fig.savefig(
    output_path,
    dpi=300,
    bbox_inches="tight",
)

plt.show()
