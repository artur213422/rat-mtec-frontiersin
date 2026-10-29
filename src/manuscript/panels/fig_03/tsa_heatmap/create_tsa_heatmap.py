"""TSA tissue GSEA heatmap."""

from pathlib import Path

import polars as pl
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from src.common.paths import PROJECT_ROOT
from src.manuscript.figures.layout.panel_layout import PanelLayout
from src.manuscript.shared.config.dataset_labels import DATASET_LABELS
from src.manuscript.shared.config.image_config import ManuscriptColors, DPI
from src.manuscript.shared.fig.fig import create_figure, save_svg, add_panel_label
from src.manuscript.shared.heatmap.draw_heatmap import draw_manuscript_heatmap
from src.manuscript.shared.heatmap.heatmap_input import ManuscriptHeatmapInput, RowOrder

INPUT_FILE: Path = PROJECT_ROOT / "data/mtec_high/derived_data/tsa_expression/samsom_et_al_tsa_analysis/gsea_per_tissue/gsea_per_tissue.csv"

NES_LIMIT: float = 2.0
P_THRESHOLD: float = 0.05
FDR_THRESHOLD: float = 0.05

def _prepare_tsa_data() -> pl.DataFrame:
    tsa_gsea: pl.DataFrame = pl.read_csv(INPUT_FILE, null_values=["NA", ""],)

    required: set[str] = {
        "dataset_id",
        "term",
        "nes",
        "nominal_p_value",
        "fdr_q_value",
    }

    missing_columns: set[str] = required - set(tsa_gsea.columns)

    if missing_columns:
        raise ValueError(f"Missing TSA GSEA columns: {sorted(missing_columns)}")

    data: pl.DataFrame = (
        tsa_gsea
        .filter(pl.col("dataset_id").is_in(list(DATASET_LABELS)))
        .select(
            pl.col("dataset_id").cast(pl.String),
            pl.col("term").cast(pl.String),
            pl.col("nes").cast(pl.Float64, strict=False),
            pl.col("nominal_p_value").cast(pl.Float64, strict=False),
            pl.col("fdr_q_value").cast(pl.Float64, strict=False),
        )
    )

    if data.is_empty():
        raise ValueError("No TSA GSEA data found for configured datasets.")

    duplicates: pl.DataFrame = data.group_by("term", "dataset_id").len().filter(pl.col("len") > 1)

    if not duplicates.is_empty():
        raise ValueError(f"Duplicate tissue/dataset results:\n{duplicates}")

    return data


def _draw_tsa_heatmap(ax: Axes, data: pl.DataFrame) -> None:
    available_datasets: set[str] = set(data.get_column("dataset_id").to_list())

    dataset_ids: list[str] = [
        dataset_id
        for dataset_id in DATASET_LABELS
        if dataset_id in available_datasets
    ]

    manuscript_input: ManuscriptHeatmapInput = ManuscriptHeatmapInput(
        data=data,
        value_col="nes",
        row_label_col="term",
        dataset_id_col="dataset_id",
        p_col="nominal_p_value",
        fdr_col="fdr_q_value",

        col_order=dataset_ids,
        col_labels=DATASET_LABELS,
        row_order=RowOrder.VALUE_ASCENDING,

        italic_row_labels=False,
        show_values=False,

        vmin=-NES_LIMIT,
        vmax=NES_LIMIT,
        min_color=ManuscriptColors.BLUE,
        max_color=ManuscriptColors.RED,
        na_color=ManuscriptColors.BLACK,

        p_threshold=P_THRESHOLD,
        fdr_threshold=FDR_THRESHOLD,

        colorbar_label="NES",
    )

    draw_manuscript_heatmap(ax, manuscript_input)


def create_tsa_heatmap(layout: PanelLayout, output_file: Path) -> Path:
    fig: Figure = create_figure(
        width_mm=layout.width,
        height_mm=layout.height,
        dpi=DPI,
    )

    fig.set_layout_engine("constrained")
    ax: Axes = fig.subplots()

    data: pl.DataFrame = _prepare_tsa_data()
    _draw_tsa_heatmap(ax, data)

    for label in ax.get_xticklabels():
        label.set_rotation(60)
        label.set_ha("left")
        label.set_rotation_mode("anchor")

    add_panel_label(ax, layout.name)
    ax.texts[-1].set_position((-0.35, 1.02))

    save_svg(fig, output_file, tight=False)
    return output_file