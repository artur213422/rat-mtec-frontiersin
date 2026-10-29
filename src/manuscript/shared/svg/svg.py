"""Combine panel SVG content into a portable, editable master SVG."""
from __future__ import annotations

from pathlib import Path
import re

from lxml import etree

from src.manuscript.figures.layout.panel_layout import FigureLayout

SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"


def _namespace_ids(root, prefix: str) -> None:
    # Matplotlib repeats IDs (figure_1, axes_1, glyphs, etc.) across files.
    ids = {element.get("id"): prefix + element.get("id")
           for element in root.iter() if element.get("id")}

    def rewrite(value: str) -> str:
        return re.sub(
            r"url\(\s*(['\"]?)#([^\s)'\"]+)\1\s*\)",
            lambda match: f"url(#{ids.get(match[2], match[2])})",
            value,
        )

    for element in root.iter():
        for key, value in list(element.attrib.items()):
            if key == "id":
                element.set(key, ids[value])
            elif key in ("href", f"{{{XLINK_NS}}}href") and value.startswith("#"):
                element.set(key, "#" + ids.get(value[1:], value[1:]))
            else:
                element.set(key, rewrite(value))
        if element.tag == f"{{{SVG_NS}}}style" and element.text:
            element.text = rewrite(element.text)


def compose_svg(layout: FigureLayout, panel_files: dict[str, Path], output_file: Path) -> None:
    missing = set(layout.panels) - panel_files.keys()
    if missing:
        raise ValueError(f"Missing SVG panels: {sorted(missing)}")

    svg = etree.Element(
        f"{{{SVG_NS}}}svg",
        nsmap={None: SVG_NS, "xlink": XLINK_NS},
        width=f"{layout.canvas['width_mm']}mm",
        height=f"{layout.canvas['height_mm']}mm",
        viewBox=f"0 0 {layout.canvas['width_mm']} {layout.canvas['height_mm']}",
    )

    for index, panel_name in enumerate(layout.panels):
        panel = layout.get_panel(panel_name)
        panel_file = Path(panel_files[panel_name])
        root = etree.parse(str(panel_file), etree.XMLParser(
            resolve_entities=False, no_network=True,
        )).getroot()
        if root.tag != f"{{{SVG_NS}}}svg":
            raise ValueError(f"Panel {panel_name} is not an SVG: {panel_file}")
        if not root.get("viewBox"):
            raise ValueError(f"Panel {panel_name} needs a viewBox: {panel_file}")
        _namespace_ids(root, f"panel_{index}_")
        # A nested SVG preserves the source coordinate system, definitions and
        # embedded images while fitting it inside the layout's millimetre box.
        root.set("x", str(panel.x))
        root.set("y", str(panel.y))
        root.set("width", str(panel.width))
        root.set("height", str(panel.height))
        root.set("preserveAspectRatio", "xMidYMid meet")
        svg.append(root)

    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_bytes(etree.tostring(
        svg, encoding="utf-8", xml_declaration=True, pretty_print=True,
    ))
