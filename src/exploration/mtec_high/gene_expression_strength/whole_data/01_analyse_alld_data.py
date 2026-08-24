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

output_folder: Path = resolve_dir_from_root("data/analysis/gene_expression_strength/whole_data/with_orthologs_only")

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

for dataset_file in source_files:


    output_file: Path = output_folder / f"{dataset_file.stem}.png"

    tec_df: pl.DataFrame = pl.read_csv(
        dataset_file,
        schema_overrides={
            "ortholog_id": pl.String,
            },
        )
    all_genes: list[str] = (tec_df
        .get_column("gene_symbol").to_list()
     )

    prepared_data: PreparedOrderedData = prepare_data(all_genes, tec_df, True)

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
        figsize=(40.0, 14.0),
    )

    fig.suptitle("all genes")

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
