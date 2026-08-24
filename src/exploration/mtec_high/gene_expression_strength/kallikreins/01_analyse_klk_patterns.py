import polars as pl
from pathlib import Path

from matplotlib import pyplot as plt
from openpyxl.styles.builtins import styles

from src.common.feratures.math.fold_change import calculate_fold_changes
from src.common.feratures.math.calc_percentile import calculate_percentile
from src.common.feratures.visualization.scatter_plot.scatter_plot import create_categorical_scatter_plot
from src.common.feratures.visualization.scatter_plot.scatter_plot_types import CategoricalScatterPlotSpec, \
    CategoricalScatterStyle, CategoricalScatterData
from src.common.paths import resolve_dir_from_root, find_file_in_project
from src.exploration.gene_expression_strength.helpers.helpers import prepare_data, PreparedOrderedData

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

source_folder: Path = resolve_dir_from_root("data/analysis/gene_expression_strength/data")
source_files: list[Path] = list(source_folder.glob("*.csv"))

chemokines_file: Path = find_file_in_project("kallikreins.csv")


chemokines_df: pl.DataFrame = pl.read_csv(
    chemokines_file,
    schema_overrides={
        "ortholog_id": pl.String,
    },
)

output_folder: Path = resolve_dir_from_root("data/analysis/gene_expression_strength/per_batch/kallikreins")

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

for dataset_file in source_files:

    chemokines: list[str] = (
        chemokines_df.filter(
            (pl.col("dataset_id") == dataset_file.stem))
        .get_column("gene_symbol")
        .to_list()
    )

    output_file: Path = output_folder / f"{dataset_file.stem}.png"

    tec_df: pl.DataFrame = pl.read_csv(
        dataset_file,
        schema_overrides={
            "ortholog_id": pl.String,
            },
        )

    prepared_data: PreparedOrderedData = prepare_data(chemokines, tec_df, False)

    style: CategoricalScatterStyle = CategoricalScatterStyle()

    data: CategoricalScatterData = CategoricalScatterData(
            categories=prepared_data.categories,
            y_values=prepared_data.y_values
    )

    #print(f"{prepared_data.categories} and tissue {prepared_data.y_values}")
    spec: CategoricalScatterPlotSpec = CategoricalScatterPlotSpec(
        style=style,
        data=data
    )

    fig, ax = plt.subplots(
        figsize=(10.0, 14.0),
    )

    fig.suptitle("chemokines")

    create_categorical_scatter_plot(ax, spec)

    ##########################################################################################
    #                                       SAVE
    ##########################################################################################

    fig.tight_layout()

    fig.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight",
    )

    print(f"Saved {output_file}")
    plt.close(fig)
