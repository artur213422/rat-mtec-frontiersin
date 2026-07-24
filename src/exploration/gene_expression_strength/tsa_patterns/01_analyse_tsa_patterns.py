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

tsa_file: Path = find_file_in_project("analysed_tsa.csv")
tsa_df: pl.DataFrame = pl.read_csv(tsa_file)

tissues: list[str] = (tsa_df.get_column("tissue").drop_nulls().unique().sort().to_list())

output_folder: Path = resolve_dir_from_root("data/analysis/gene_expression_strength/per_batch/tsa_patterns_result")

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

for tissue in tissues:
    print("Running:", tissue)
    tissue_folder: Path = output_folder / tissue
    tissue_folder.mkdir(parents=True, exist_ok=True)

    tsa_genes: list[str] = (
        tsa_df.filter(
            (pl.col("tissue") == tissue) & (pl.col("found_in_dataset") == True))
        .get_column("gene_symbol")
        .drop_nulls()
        .unique(maintain_order=True)
        .to_list()
    )

    for dataset_file in source_files:

        output_file: Path = tissue_folder / f"{dataset_file.stem}.png"

        tec_df: pl.DataFrame = pl.read_csv(
            dataset_file,
            schema_overrides={
                "ortholog_id": pl.String,
            },
        )

        prepared_data: PreparedOrderedData = prepare_data(tsa_genes, tec_df, True)

        if len(prepared_data.categories) == 0:
            print(f"Skipping {dataset_file.stem} and tissue {tissue} because there are no genes")
            continue

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
        fig.suptitle(tissue)

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
