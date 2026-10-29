from pathlib import Path

import cairosvg


def save_png(svg_file: Path, png_file: Path) -> None:
    png_file.parent.mkdir(parents=True, exist_ok=True)

    cairosvg.svg2png(
        bytestring=svg_file.read_bytes(),
        write_to=str(png_file),
        dpi=300,
    )