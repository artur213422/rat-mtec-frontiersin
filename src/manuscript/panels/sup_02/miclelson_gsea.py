from pathlib import Path

from src.manuscript.figures.layout.panel_layout import PanelLayout
from src.manuscript.panels.sup_02.gsea_heatmap import SOURCES, create_gsea_heatmap


def create_michelson_gsea(layout: PanelLayout, output_file: Path) -> Path:
    return create_gsea_heatmap(layout, output_file, SOURCES[0])
