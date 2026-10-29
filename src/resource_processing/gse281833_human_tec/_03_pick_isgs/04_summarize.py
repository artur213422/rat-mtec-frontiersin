from pathlib import Path

import polars as pl

from src.common.paths import resolve_dir_from_root

picked_isgs_dir: Path = resolve_dir_from_root(
    "resources/interferon/gse28183_tec/standardized/picked_isgs"
)
summary_file: Path = picked_isgs_dir / "picked_isgs_summary.txt"

summary_lines: list[str] = []

for isgs_file in sorted(picked_isgs_dir.glob("*.csv")):
    isgs_df: pl.DataFrame = pl.read_csv(isgs_file)

    summary_lines.extend([
        isgs_file.name,
        f"Total candidates: {isgs_df.height}",
    ])

    counts_per_ifn = (
        isgs_df
        .group_by("target_condition")
        .len()
        .sort("target_condition")
    )

    for row in counts_per_ifn.iter_rows(named=True):
        summary_lines.append(
            f"  {row['target_condition']}: {row['len']}"
        )

    summary_lines.append("")

summary_file.write_text("\n".join(summary_lines), encoding="utf-8")