from pathlib import Path

import numpy as np
import polars as pl
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from src.common.dataset_catalog.api.get_datasets import read_file_with_ortholog_id
from src.common.paths import PROJECT_ROOT
from src.manuscript.figures.layout.panel_layout import PanelLayout
from src.manuscript.shared.config.dataset_labels import DATASET_LABELS
from src.manuscript.shared.config.image_config import ManuscriptColors, DPI
from src.manuscript.shared.fig.fig import create_figure, save_svg, add_panel_label
from src.manuscript.shared.heatmap.draw_heatmap import draw_manuscript_heatmap
from src.manuscript.shared.heatmap.heatmap_input import ManuscriptHeatmapInput, RowOrder

INPUT_FILE: Path = PROJECT_ROOT / "data/mtec_high/derived_data/chemokine_expression/chemokines.csv"

P_THRESHOLD: float = 0.05
FDR_THRESHOLD: float = 0.05

def _prepare_chemokine_data() -> pl.DataFrame:
    chemokines: pl.DataFrame = read_file_with_ortholog_id(INPUT_FILE)

    required: set[str] = {
        "gene_symbol",
        "dataset_id",
        "log2_fold_change",
        "p_val",
        "p_adj",
    }

    missing_columns: set[str] = required - set(chemokines.columns)

    if missing_columns:
        raise ValueError(f"Missing chemokine columns: {sorted(missing_columns)}")

    data: pl.DataFrame = (
        chemokines
        .filter(pl.col("dataset_id").is_in(list(DATASET_LABELS)))
        .select(
            pl.col("dataset_id").cast(pl.String),

            pl.when(pl.col("dataset_id").str.starts_with("mouse_"))
            .then(pl.col("gene_symbol").str.replace(r"_family$", ""))
            .otherwise(pl.col("gene_symbol"))
            .cast(pl.String)
            .alias("gene_symbol"),

            pl.col("log2_fold_change").cast(pl.Float64),
            pl.col("p_val").cast(pl.Float64),
            pl.col("p_adj").cast(pl.Float64),
        )
    )

    if data.is_empty():
        raise ValueError("No chemokine data found for configured datasets.")

    duplicates: pl.DataFrame = data.group_by("gene_symbol", "dataset_id").len().filter(pl.col("len") > 1)

    if not duplicates.is_empty():
        raise ValueError(f"Duplicate gene/dataset results:\n{duplicates}")

    return data


def _draw_chemokine_heatmap(ax: Axes, data: pl.DataFrame) -> None:
    available_datasets: set[str] = set(data.get_column("dataset_id").to_list())

    dataset_ids: list[str] = [
        dataset_id
        for dataset_id in DATASET_LABELS
        if dataset_id in available_datasets
    ]

    values: np.ndarray = data.get_column("log2_fold_change").to_numpy()
    finite_values: np.ndarray = values[np.isfinite(values)]

    if finite_values.size == 0:
        raise ValueError("Chemokine heatmap contains no finite fold changes.")

    vmin: float = min(float(finite_values.min()), -0.01)
    vmax: float = max(float(finite_values.max()), 0.01)

    manuscript_input: ManuscriptHeatmapInput = ManuscriptHeatmapInput(
        data=data,
        value_col="log2_fold_change",
        row_label_col="gene_symbol",
        dataset_id_col="dataset_id",
        p_col="p_val",
        fdr_col="p_adj",

        col_order=dataset_ids,
        col_labels=DATASET_LABELS,
        row_order=RowOrder.VALUE_ASCENDING,

        italic_row_labels=True,
        show_values=False,

        vmin=vmin,
        vmax=vmax,
        min_color=ManuscriptColors.GREEN,
        max_color=ManuscriptColors.RED,
        na_color=ManuscriptColors.BLACK,

        p_threshold=P_THRESHOLD,
        fdr_threshold=FDR_THRESHOLD,

        colorbar_label="log₂ fold change (KO/control)",
    )

    draw_manuscript_heatmap(ax, manuscript_input)

def _recover_chemokines(data: pl.DataFrame) -> pl.DataFrame:
    return data


def create_chemokine_panel(layout: PanelLayout, output_file: Path) -> Path:
    fig: Figure = create_figure(
        width_mm=layout.width,
        height_mm=layout.height,
        dpi=DPI,
    )

    fig.set_layout_engine("constrained")
    ax: Axes = fig.subplots()

    data: pl.DataFrame = _prepare_chemokine_data()
    data = _recover_chemokines(data)
    _draw_chemokine_heatmap(ax, data)

    for label in ax.get_xticklabels():
        label.set_rotation(60)
        label.set_ha("left")
        label.set_rotation_mode("anchor")

    add_panel_label(ax, layout.name)
    ax.texts[-1].set_position((-0.35, 1.02))

    save_svg(fig, output_file, tight=False)
    return output_file