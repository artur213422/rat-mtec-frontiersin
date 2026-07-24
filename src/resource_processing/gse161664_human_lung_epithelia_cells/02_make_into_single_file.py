import re
from pathlib import Path

import polars as pl

from src.common.paths import resolve_dir_from_root


def get_sample_name(file_path: Path) -> str:
    match = re.search(
        pattern=r"_SAEC_(?P<treatment>[^_]+)_rep(?P<rep>\d+)$",
        string=file_path.stem,
    )

    if match is None:
        raise ValueError(f"Cannot parse sample name: {file_path.name}")

    treatment: str = match.group("treatment")
    replicate: str = match.group("rep")

    treatment_names: dict[str, str] = {
        "Control": "ctrl",
        "IFNa": "ifna",
        "IFNb": "ifnb",
        "IFNg": "ifng",
        "IFNl3": "ifnl3",
    }

    if treatment not in treatment_names:
        raise ValueError(f"Unknown treatment: {treatment}")

    return f"{treatment_names[treatment]}{replicate}"


output_dir: Path = resolve_dir_from_root("resources/interferon/gse161664/raw")

paths: list[Path] = sorted(output_dir.glob("*.txt"))

dfs: list[pl.DataFrame] = []

for file_path in paths:
    sample_name: str = get_sample_name(file_path)

    sample_df: pl.DataFrame = (
        pl.read_csv(
            file_path,
            separator="\t",
            has_header=False,
            new_columns=["gene", sample_name],
        )
        .filter(
            ~pl.col("gene").str.starts_with("__")
        )
    )

    dfs.append(sample_df)

result_df: pl.DataFrame = (
    pl.concat(
        dfs,
        how="align_full",
    )
    .sort("gene")
)

result_path: Path = resolve_dir_from_root("resources/interferon/gse161664/ready_to_use") / "gse161664_counts.csv"
result_df.write_csv(result_path)

print(result_df)
print(f"Saved to: {result_path}")