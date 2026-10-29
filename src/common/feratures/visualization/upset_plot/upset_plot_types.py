from dataclasses import dataclass, field
from typing import TypeAlias

from src.common.types.enums import Colors


CategoryName: TypeAlias = str


@dataclass(frozen=True, slots=True)
class UpSetBarSegment:
    count: int
    color: Colors
    label: str | None = None

    def __post_init__(self) -> None:
        if self.count < 0:
            raise ValueError(f"Segment count cannot be negative: {self.count}")


@dataclass(frozen=True, slots=True)
class UpSetIntersection:
    tuple_label: str
    segments: tuple[UpSetBarSegment, ...]

    def __post_init__(self) -> None:
        if not self.segments:
            raise ValueError("An intersection must contain at least one bar segment.")

    @property
    def total_count(self) -> int:
        return sum(
            segment.count
            for segment in self.segments
        )


@dataclass(frozen=True, slots=True)
class UpSetPlotData:
    intersections: tuple[UpSetIntersection, ...]

@dataclass(frozen=True, slots=True)
class UpSetPlotStyle:
    bar_width: float = 0.7
    dot_size: float = 50.0
    inactive_dot_size: float = 25.0
    line_width: float = 2.0
    show_bar_counts: bool = True
    matrix_height_fraction: float = 0.35
    x_tick_rotation: float = 45.0

    collapse_duplicate_tuple_labels: bool = False

@dataclass(frozen=True, slots=True)
class UpSetPlotSpec:
    data: UpSetPlotData
    style: UpSetPlotStyle = field(
        default_factory=UpSetPlotStyle
    )