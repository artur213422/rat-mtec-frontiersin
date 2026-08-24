from pathlib import Path

import polars as pl

def summarize(data: pl.DataFrame, summary_file: Path):
    summary_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    ifns: list[str] = (
        data
        .get_column("ifn")
        .unique()
        .sort()
        .to_list()
    )

    summary_lines: list[str] = []

    for ifn in ifns:
        ifn_data: pl.DataFrame = data.filter(
            pl.col("ifn") == ifn
        )

        total_count: int = ifn_data.height

        unique_count: int = (
            ifn_data
            .filter(pl.col("has_clean_ortholog"))
            .height
        )

        summary_lines.extend(
            [
                f"{ifn} summary:",
                f"Unique: {unique_count}",
                f"Total: {total_count}",
                "",
            ]
        )

    summary_file.write_text(
        "\n".join(summary_lines),
        encoding="utf-8",
    )

    print(f"Summary saved to: {summary_file}")