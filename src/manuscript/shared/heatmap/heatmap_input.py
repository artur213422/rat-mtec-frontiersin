from dataclasses import dataclass
from enum import StrEnum
import polars as pl

from src.manuscript.shared.config.image_config import ManuscriptColors

class RowOrder(StrEnum):
    ALPHABETICAL = "alphabetical"
    VALUE_ASCENDING = "value_ascending"

@dataclass(frozen=True, slots=True)
class ManuscriptHeatmapInput:
    # Data
    data: pl.DataFrame
    value_col: str      # NES, fold change
    row_label_col: str  # gene_symbol, pathway
    dataset_id_col: str
    p_col: str
    fdr_col: str

    # Ordering
    col_order: list[str]
    row_order: RowOrder | list[str]

    # Visual
    italic_row_labels: bool = False
    show_values: bool = False
    col_labels: dict[str, str] | None = None

    # Colors
    vmin: float = -2.0
    vmax: float = 2.0
    min_color: str = ManuscriptColors.BLUE
    max_color: str = ManuscriptColors.RED
    na_color: str = ManuscriptColors.WHITE

    # Legend
    colorbar_label: str = ""

    # Significance
    p_threshold: float = 0.05
    fdr_threshold: float = 0.05

