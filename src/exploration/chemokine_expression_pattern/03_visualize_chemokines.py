from pathlib import Path
import polars as pl
from matplotlib import pyplot as plt

from src.common.feratures.visualization.heatmap.heatmap import create_heatmap
from src.common.feratures.visualization.heatmap.heatmap_models import HeatmapStyle, HeatmapSpec, HeatmapData
from src.common.models.dataset_metadata.model import DatasetCatalog
from src.common.models.dataset_metadata.ops import build_dataset_catalog
from src.common.paths import try_find_file_in_project, resolve_dir_from_root
from src.common.types.enums import Colors

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


output_file: Path = resolve_dir_from_root("data/analysis/chemokine_expression_pattern/results") / "chemokines.png"

analysed_file: Path = try_find_file_in_project("analysed_chemokines.csv")

analysed_tsa_df: pl.DataFrame = pl.read_csv(
    analysed_file,
    schema_overrides={
        "ortholog_id": pl.String,
    },
)

datasets: DatasetCatalog = build_dataset_catalog()

rat_dataset: str = datasets.get_rat_dataset().dataset_id

##########################################################################################
#                                       PROCESS DATA
##########################################################################################




heatmap_data: HeatmapData = HeatmapData.from_long_frame(
    data=analysed_tsa_df,
    row_label_col="gene_symbol",
    col_label_col="dataset_id",
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
        # vmax=3,
        # vmin=-7,

 )

spec: HeatmapSpec = HeatmapSpec(
        data=heatmap_data,
        style=style
)

fig, ax = plt.subplots(
    figsize=(10, 14),
)
fig.suptitle(f"Chemokines")

ax = create_heatmap(ax, spec)

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


