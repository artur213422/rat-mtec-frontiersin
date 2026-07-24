import math
from dataclasses import dataclass
from src.common.types.enums import Colors


@dataclass(frozen=True, slots=True)
class CategoricalScatterStyle:
    marker_color: Colors = Colors.BLUE
    marker_size: float = 20.0
    marker_alpha: float = 0.8


@dataclass(slots=True)
class CategoricalScatterData:
    categories: list[str]
    y_values: list[float]

    def __post_init__(self) -> None:
        if not self.categories:
            raise ValueError("Categories cannot be empty.")

        if len(self.categories) != len(self.y_values):
            raise ValueError("Categories and y_values must have the same length. Received {len(self.categories)} categories and {len(self.y_values)} y-values.")

        if any(not math.isfinite(value) for value in self.y_values):
            raise ValueError("All y_values must be finite numeric values.")

@dataclass(frozen=True, slots=True)
class CategoricalScatterPlotSpec:
    style: CategoricalScatterStyle
    data: CategoricalScatterData

