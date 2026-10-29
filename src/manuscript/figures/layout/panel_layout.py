from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json


@dataclass(frozen=True)
class PanelLayout:
    name: str
    x: float
    y: float
    width: float
    height: float


class FigureLayout:

    def __init__(self, json_file: Path):

        data = json.loads(
            json_file.read_text()
        )

        self.canvas = data["canvas"]

        self.panels: dict[str, PanelLayout] = {}

        for name, panel in data["panels"].items():

            self.panels[name] = PanelLayout(
                name=name,
                x=panel["x"],
                y=panel["y"],
                width=panel["width"],
                height=panel["height"],
            )


    def get_panel(self, name: str) -> PanelLayout:

        return self.panels[name]