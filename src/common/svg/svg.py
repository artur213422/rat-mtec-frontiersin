from pathlib import Path

from src.common.paths import PROJECT_ROOT
import cairosvg

source = PROJECT_ROOT / "data/mtec_high/derived_data/standard_gsea/simplify_enrichment"
files: list[Path] = sorted(source.rglob("*.svg"))

print(len(files))

for svg_file in files:
    png_file = svg_file.with_suffix(".png")
    print(f"Converting {svg_file}")

    cairosvg.svg2png(
        bytestring=svg_file.read_bytes(),
        url=svg_file.as_uri(),
        write_to=str(png_file),
    )