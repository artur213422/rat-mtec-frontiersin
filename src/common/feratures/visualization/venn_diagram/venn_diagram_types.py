from dataclasses import dataclass, field

from src.common.types.enums import Colors


@dataclass(frozen=True, slots=True)
class VennDiagramSet:
    label: str
    data: frozenset[str]
    color: Colors

    def __post_init__(self) -> None:
        if not self.label.strip():
            raise ValueError("Venn diagram set label cannot be empty.")

        normalized_data: frozenset[str] = frozenset(
            value.strip().lower()
            for value in self.data
            if value.strip()
        )

        object.__setattr__(
            self,
            "data",
            normalized_data,
        )

@dataclass(frozen=True, slots=True)
class VennDiagramData:
    sets: tuple[VennDiagramSet, ...]

    def __post_init__(self) -> None:
        if len(self.sets) not in (2, 3):
            raise ValueError("A Venn diagram requires two or three sets.")

@dataclass(frozen=True, slots=True)
class VennDiagramStyle:
    alpha: float = 0.5


@dataclass(frozen=True, slots=True)
class VennDiagramSpec:
    data: VennDiagramData
    style: VennDiagramStyle = field(
        default_factory=VennDiagramStyle
    )