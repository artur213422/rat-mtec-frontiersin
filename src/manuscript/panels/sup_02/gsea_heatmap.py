"""Shared mTEC high cell-signature GSEA panels for Supplementary Figure 2."""

from pathlib import Path

import polars as pl

from src.common.paths import PROJECT_ROOT
from src.manuscript.figures.layout.panel_layout import PanelLayout
from src.manuscript.shared.config.dataset_labels import DATASET_LABELS
from src.manuscript.shared.config.image_config import DPI, TITLE_SIZE, ManuscriptColors
from src.manuscript.shared.fig.fig import create_figure, save_svg
from src.manuscript.shared.heatmap.draw_heatmap import draw_manuscript_heatmap
from src.manuscript.shared.heatmap.heatmap_input import ManuscriptHeatmapInput, RowOrder

SOURCES = [
    ("Michelson et al., 2022", "michelson_2022_mimetic_cells", "gsea_analysis.csv"),
    ("Baran-Gale et al., 2020", "hollander_thymus_cell_types_gsea", "cell_type_gsea.csv"),
]
IGNORED_TERMS = {
    "Baran-Gale et al., 2020": ["Mature cTEC", "Perinatal cTEC"],
}
NES_LIMIT = 3.5


def load_gsea_data(source: tuple[str, str, str]) -> pl.DataFrame:
    title, folder, filename = source
    path = PROJECT_ROOT / "data/mtec_high/derived_data" / folder / filename
    data = pl.read_csv(path, null_values=["NA", ""])
    required = {"dataset_id", "term", "nes", "nominal_p_value", "fdr_q_value"}
    missing = required - set(data.columns)
    if missing:
        raise ValueError(f"Missing GSEA columns in {path}: {sorted(missing)}")
    data = data.filter(
        pl.col("dataset_id").str.ends_with("_mtec_high"),
        ~pl.col("term").is_in(IGNORED_TERMS.get(title, [])),
    ).select(
        pl.col("dataset_id").cast(pl.String),
        pl.col("term").cast(pl.String),
        pl.col("nes").cast(pl.Float64),
        pl.col("nominal_p_value").cast(pl.Float64),
        pl.col("fdr_q_value").cast(pl.Float64),
    )
    missing_datasets = set(DATASET_LABELS) - set(data["dataset_id"].to_list())
    if missing_datasets:
        raise ValueError(f"Missing mTEC high datasets in {path}: {sorted(missing_datasets)}")
    return data


def create_gsea_heatmap(
    layout: PanelLayout, output_file: Path, source: tuple[str, str, str],
) -> Path:
    data = load_gsea_data(source)
    available = set(data["dataset_id"].to_list())
    dataset_ids = list(DATASET_LABELS) + sorted(available - set(DATASET_LABELS))
    fig = create_figure(layout.width, layout.height, dpi=DPI)
    # Reserve space for row labels, rotated dataset labels, and inset legends.
    ax = fig.add_axes((0.30, 0.12, 0.52, 0.63))
    draw_manuscript_heatmap(
        ax,
        ManuscriptHeatmapInput(
            data=data,
            value_col="nes",
            row_label_col="term",
            dataset_id_col="dataset_id",
            p_col="nominal_p_value",
            fdr_col="fdr_q_value",
            col_order=dataset_ids,
            col_labels=DATASET_LABELS,
            row_order=RowOrder.VALUE_ASCENDING,
            show_values=False,
            vmin=-NES_LIMIT,
            vmax=NES_LIMIT,
            na_color=ManuscriptColors.BLACK,
            colorbar_label="NES",
            p_threshold=0.05,
            fdr_threshold=0.05,
        ),
    )
    for label in ax.get_xticklabels():
        label.set_rotation(60)
        label.set_ha("left")
        label.set_rotation_mode("anchor")
    fig.text(0.03, 0.96, layout.name, fontsize=12, fontweight="bold", va="top")
    fig.text(0.30, 0.96, source[0], fontsize=TITLE_SIZE, va="top")
    return save_svg(fig, output_file, tight=False)
