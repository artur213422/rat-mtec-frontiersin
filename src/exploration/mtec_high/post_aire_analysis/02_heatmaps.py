from pathlib import Path

import polars as pl
from matplotlib import pyplot as plt

from src.common.dataset_catalog.api.get_datasets import get_mtec_low, get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.feratures.visualization.heatmap.heatmap import create_heatmap
from src.common.feratures.visualization.heatmap.heatmap_models import HeatmapData, HeatmapSpec, HeatmapStyle

from src.common.paths import find_file_in_project, resolve_dir_from_root
from src.common.types.enums import Colors


##########################################################################################
#                                       SETTINGS
##########################################################################################

MAX_GENES_PER_CELL_TYPE: int = 400

MARKER_SUBTYPE_COL: str = "Subtype"
MARKER_GENE_COL: str = "GeneName"
MARKER_P_VALUE_COL: str = "PV(adj)"

DATASET_GENE_COL: str = "gene_symbol"
DATASET_VALUE_COL: str = "log2_fold_change"


##########################################################################################
#                                       PREPARE DATA
##########################################################################################

marker_file: Path = find_file_in_project("thymus_cell_markers_hollander_2020.csv")

output_folder: Path = resolve_dir_from_root("data/mtec_high/analysis/post_aire_analysis")

output_folder.mkdir(parents=True, exist_ok=True)

datasets: list[DatasetMetadata] = get_mtec_high()


##########################################################################################
#                               PREPARE CELL-TYPE MARKERS
##########################################################################################

marker_df: pl.DataFrame = (
    pl.read_csv(
        marker_file,
        separator=";",
    )
    .select(
        pl.col(MARKER_SUBTYPE_COL).alias("cell_type"),

        # Normalize case because the standardized expression tables
        # use lower-case gene symbols.
        pl.col(MARKER_GENE_COL)
        .str.to_lowercase()
        .alias("gene_symbol"),

        pl.col(MARKER_P_VALUE_COL)
        .cast(pl.Float64)
        .alias("marker_p_adj"),
    )
    .filter(
        pl.col("gene_symbol").is_not_null()
        & (pl.col("gene_symbol") != "")
        & (pl.col("gene_symbol") != "na")
    )
    .unique(
        subset=["cell_type", "gene_symbol"],
        keep="first",
    )
    .sort(
        ["cell_type", "marker_p_adj"]
    )
    .with_columns(
        pl.col("marker_p_adj")
        .rank(method="ordinal")
        .over("cell_type")
        .alias("marker_rank")
    )
    .filter(
        pl.col("marker_rank") <= MAX_GENES_PER_CELL_TYPE
    )
)


print(
    marker_df
    .group_by("cell_type")
    .agg(
        pl.len().alias("selected_genes"),
        pl.col("marker_p_adj").max().alias("worst_selected_p_adj"),
    )
    .sort("cell_type")
)


##########################################################################################
#                           COLLECT ALL DATASETS INTO LONG TABLE
##########################################################################################

dataset_frames: list[pl.DataFrame] = []


for dataset in datasets:

    print(f"Processing dataset: {dataset.dataset_id}")

    tec_df: pl.DataFrame = dataset.get_standardized_df()

    selected_df: pl.DataFrame = (
        tec_df
        .select(
            pl.col(DATASET_GENE_COL)
            .str.to_lowercase()
            .alias("gene_symbol"),

            pl.col(DATASET_VALUE_COL),
        )
        .drop_nulls(
            ["gene_symbol", DATASET_VALUE_COL]
        )
        .with_columns(
            pl.lit(dataset.dataset_id).alias("dataset")
        )
    )

    # HeatmapData requires one value per gene × dataset.
    #
    # This should normally already be true in standardized datasets.
    # If multiple probes/rows remain for a gene, mean FC prevents
    # duplicate pivot cells.
    selected_df = (
        selected_df
        .group_by(
            ["gene_symbol", "dataset"],
            maintain_order=True,
        )
        .agg(
            pl.col(DATASET_VALUE_COL).mean()
        )
    )

    dataset_frames.append(selected_df)


all_expression_df: pl.DataFrame = pl.concat(
    dataset_frames,
    how="vertical",
)


##########################################################################################
#                     JOIN CELL-TYPE MARKERS TO EXPRESSION DATA
##########################################################################################

heatmap_source_df: pl.DataFrame = (
    marker_df
    .join(
        all_expression_df,
        on="gene_symbol",
        how="inner",
    )
)


##########################################################################################
#                               PRINT COVERAGE DIAGNOSTICS
##########################################################################################

coverage_df: pl.DataFrame = (
    heatmap_source_df
    .group_by(
        ["cell_type", "dataset"]
    )
    .agg(
        pl.col("gene_symbol")
        .n_unique()
        .alias("genes_found")
    )
    .sort(
        ["cell_type", "dataset"]
    )
)

print("\nMarker coverage:")
print(coverage_df)


##########################################################################################
#                                   CREATE HEATMAPS
##########################################################################################

cell_types: list[str] = (
    marker_df
    .get_column("cell_type")
    .unique(maintain_order=True)
    .to_list()
)


for cell_type in cell_types:

    selected_data: pl.DataFrame = (
        heatmap_source_df
        .filter(
            pl.col("cell_type") == cell_type
        )
        .select(
            "gene_symbol",
            "dataset",
            DATASET_VALUE_COL,
        )
    )

    if selected_data.is_empty():
        print(f"{cell_type}: no genes found. Skipping.")
        continue

    n_genes: int = (
        selected_data
        .get_column("gene_symbol")
        .n_unique()
    )

    print(
        f"Creating {cell_type}: "
        f"{n_genes} genes represented"
    )


    ######################################################################################
    #                                   HEATMAP DATA
    ######################################################################################

    heatmap_data: HeatmapData = HeatmapData.from_long_frame(
        data=selected_data,
        row_label_col="gene_symbol",
        col_label_col="dataset",
        value_col=DATASET_VALUE_COL,
    )


    ######################################################################################
    #                           MOVE RAT LEFT / SORT BY RAT
    ######################################################################################

    rat_dataset_ids: list[str] = [
        dataset.dataset_id
        for dataset in datasets
        if "rat" in dataset.dataset_id.lower()
    ]

    if len(rat_dataset_ids) == 1:
        rat_dataset_id: str = rat_dataset_ids[0]

        heatmap_data.move_col_to_left(rat_dataset_id)
        heatmap_data.sort_rows_by_col(rat_dataset_id)


    ######################################################################################
    #                                       STYLE
    ######################################################################################

    style: HeatmapStyle = HeatmapStyle(
        high_color=Colors.RED,
        low_color=Colors.GREEN,
        center_color=Colors.WHITE,
        nan_cell_color=Colors.BLACK,
        vcenter=0,
        # vmax=3,
        # vmin=-3,
    )

    spec: HeatmapSpec = HeatmapSpec(
        data=heatmap_data,
        style=style,
    )


    ######################################################################################
    #                                       PLOT
    ######################################################################################

    # 500 labelled genes need a tall figure.
    # Cap the height so the PNG does not become absurdly large.
    figure_height: float = min(
        50.0,
        max(
            8.0,
            n_genes * 0.10,
        ),
    )

    fig, ax = plt.subplots(
        figsize=(8, figure_height)
    )

    fig.suptitle(
        f"{cell_type} markers"
    )

    create_heatmap(
        ax,
        spec,
    )


    ######################################################################################
    #                                       SAVE
    ######################################################################################

    safe_cell_type: str = (
        cell_type
        .replace(" ", "_")
        .replace("/", "_")
    )

    output_path: Path = (
        output_folder
        / f"{safe_cell_type}.png"
    )

    fig.tight_layout()

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


print("Done")