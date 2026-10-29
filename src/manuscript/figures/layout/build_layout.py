from __future__ import annotations

import argparse
from collections import deque
import json
import math
from pathlib import Path

from PIL import Image


LAYOUT_DIR = Path(__file__).resolve().parent
fig_1_png = LAYOUT_DIR / "figure_1.png"
fig_2_png = LAYOUT_DIR / "figure_2.png"
fig_3_png = LAYOUT_DIR / "figure_3.png"
sup_1_png = LAYOUT_DIR / "sup_1.png"
sup_2_png = LAYOUT_DIR / "sup_2.png"

fig_1_json = LAYOUT_DIR / "figure_1_layout.json"
fig_2_json = LAYOUT_DIR / "figure_2_layout.json"
fig_3_json = LAYOUT_DIR / "figure_3_layout.json"
sup_1_json = LAYOUT_DIR / "sup_1_layout.json"
sup_2_json = LAYOUT_DIR / "sup_2_layout.json"


def _panel_label(index: int) -> str:
    """Return spreadsheet-style labels for a zero-based index."""
    label = ""
    index += 1
    while index:
        index, remainder = divmod(index - 1, 26)
        label = chr(ord("A") + remainder) + label
    return label


def _rectangles(image: Image.Image) -> list[tuple[int, int, int, int, str]]:
    """Find four-connected regions of equal RGB color and require rectangles."""
    width, height = image.size
    pixels = list(image.convert("RGBA").getdata())
    visited = bytearray(width * height)
    rectangles = []

    for start, pixel in enumerate(pixels):
        if visited[start]:
            continue
        visited[start] = 1
        color = pixel[:3]
        if pixel[3] == 0 or color == (255, 255, 255):
            continue

        queue = deque([start])
        left = right = start % width
        top = bottom = start // width
        area = 0
        while queue:
            position = queue.popleft()
            x, y = position % width, position // width
            left, right = min(left, x), max(right, x)
            top, bottom = min(top, y), max(bottom, y)
            area += 1
            neighbors = []
            if x:
                neighbors.append(position - 1)
            if x + 1 < width:
                neighbors.append(position + 1)
            if y:
                neighbors.append(position - width)
            if y + 1 < height:
                neighbors.append(position + width)
            for neighbor in neighbors:
                candidate = pixels[neighbor]
                if not visited[neighbor] and candidate[3] != 0 and candidate[:3] == color:
                    visited[neighbor] = 1
                    queue.append(neighbor)

        box_width, box_height = right - left + 1, bottom - top + 1
        hex_color = "#{:02X}{:02X}{:02X}".format(*color)
        if area != box_width * box_height:
            raise ValueError(
                f"Region {hex_color} at pixel ({left}, {top}) is not a solid "
                "rectangle. Use flat fills without outlines or antialiasing."
            )
        rectangles.append((left, top, box_width, box_height, hex_color))

    if not rectangles:
        raise ValueError("The PNG contains no non-white, non-transparent panels.")
    return sorted(rectangles, key=lambda box: (box[1], box[0]))


def build_layout(
    png_path: Path,
    json_path: Path | None = None,
    *,
    width_mm: float = 180.0,
    height_mm: float | None = None,
) -> dict:
    """Write a layout, preserving PNG proportions unless height is explicit.

    The default physical width is 180 mm; height follows the PNG aspect ratio.
    """
    sizes = (width_mm,) if height_mm is None else (width_mm, height_mm)
    if not all(math.isfinite(size) and size > 0 for size in sizes):
        raise ValueError("Canvas dimensions must be finite, positive millimetres.")
    png_path = Path(png_path)
    json_path = (
        Path(json_path) if json_path is not None
        else png_path.with_name(f"{png_path.stem}_layout.json")
    )
    if png_path.resolve() == json_path.resolve():
        raise ValueError("The JSON output must differ from the input PNG.")
    with Image.open(png_path) as image:
        if image.format != "PNG":
            raise ValueError(f"Expected a PNG image: {png_path}")
        pixel_width, pixel_height = image.size
        rectangles = _rectangles(image)

    if height_mm is None:
        height_mm = width_mm * pixel_height / pixel_width

    existing_panels = {}
    if json_path.exists():
        with json_path.open(encoding="utf-8") as handle:
            existing_panels = json.load(handle).get("panels", {})

    panels = {}
    for index, (left, top, width, height, color) in enumerate(rectangles):
        label = _panel_label(index)
        panels[label] = {
            "name": existing_panels.get(label, {}).get("name", label),
            "x": left * width_mm / pixel_width,
            "y": top * height_mm / pixel_height,
            "width": width * width_mm / pixel_width,
            "height": height * height_mm / pixel_height,
            "color": color,
        }
    layout = {
        "canvas": {
            "width_mm": width_mm,
            "height_mm": height_mm,
            "origin": "top-left",
        },
        "order": list(panels),
        "panels": panels,
    }
    with json_path.open("w", encoding="utf-8") as handle:
        json.dump(layout, handle, indent=2, allow_nan=False)
        handle.write("\n")
    return layout


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "pngs", type=Path, nargs="*", default=[fig_1_png, fig_2_png, fig_3_png, sup_1_png, sup_2_png],
        help="PNG layouts; defaults to figure_1.png, figure_2.png and figure_3.png beside this script.",
    )
    parser.add_argument("--width-mm", type=float, default=180.0)
    parser.add_argument(
        "--height-mm", type=float, default=None,
        help="Physical height; defaults to the PNG aspect ratio at --width-mm.",
    )
    args = parser.parse_args()
    for png_path in args.pngs:
        try:
            layout = build_layout(png_path, width_mm=args.width_mm, height_mm=args.height_mm)
        except (OSError, ValueError) as error:
            parser.exit(1, f"{png_path}: {error}\n")
        output = png_path.with_name(f"{png_path.stem}_layout.json")
        print(f"{output}: {len(layout['panels'])} panels")


if __name__ == "__main__":
    main()
