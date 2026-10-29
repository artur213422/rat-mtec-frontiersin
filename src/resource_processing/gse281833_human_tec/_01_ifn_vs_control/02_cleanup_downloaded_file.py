import re
from pathlib import Path
from typing import Final

import polars as pl
from src.common.paths import PROJECT_ROOT

##########################################################################################
#                                       PREPARE DATA
##########################################################################################
RESOURCE_DIR: Final[Path] = PROJECT_ROOT / "resources/interferon/gse28183_tec"

SOURCE_FILE: Final[Path] = RESOURCE_DIR / "source" / "GSE281833_count_matrix.tsv"


OUTPUT_FILE: Final[Path] = RESOURCE_DIR / "intermediate" / "counts" / "gse281833_counts.csv"



COLUMN_PATTERN: Final[re.Pattern[str]] = re.compile(r"^ZE_(Ctrl|IFNA|IFNG)_(\d+)\.(\d+)$")

TREATMENT_ORDER: Final[dict[str, int]] = {
    "Ctrl": 0,
    "IFNA": 1,
    "IFNG": 2,
}

TREATMENT_NAMES: Final[dict[str, str]] = {
    "Ctrl": "ctrl",
    "IFNA": "ifna",
    "IFNG": "ifng",
}


def parse_sample_column(column: str) -> tuple[str, int, int] | None:
    match: re.Match[str] | None = COLUMN_PATTERN.fullmatch(column)

    if match is None:
        return None

    treatment: str = match.group(1)
    experiment: int = int(match.group(2))
    replicate: int = int(match.group(3))

    return treatment, experiment, replicate


##########################################################################################
#                                       PROCESS DATA
##########################################################################################

def prepare_count_matrix(source_file: Path,  output_file: Path,) -> pl.DataFrame:
    dataframe: pl.DataFrame = pl.read_csv(source_file,  separator="\t")

    sample_columns: list[tuple[str, str, int, int]] = []

    for column in dataframe.columns:
        parsed: tuple[str, int, int] | None = parse_sample_column(column)

        if parsed is None:
            continue

        treatment, experiment, replicate = parsed
        sample_columns.append(
            (column, treatment, experiment, replicate)
        )

    # Order by treatment, then experiment, then replicate.
    # This places experiment 1 before experiment 2.
    sample_columns.sort(
        key=lambda item: (
            TREATMENT_ORDER[item[1]],
            item[2],
            item[3],
        )
    )

    if not sample_columns:
        raise ValueError(f"No control/IFNA/IFNG sample columns found in {source_file}")
    for treatment in TREATMENT_ORDER:
        if not any(item[1] == treatment for item in sample_columns):
            raise ValueError(f"Missing {treatment} samples in {source_file}")

    selected_columns: list[str] = [
        "gene_id",
        *(column for column, _, _, _ in sample_columns),
    ]

    rename_mapping: dict[str, str] = {
        "gene_id": "gene_symbol",
    }

    treatment_counts: dict[str, int] = {
        "Ctrl": 0,
        "IFNA": 0,
        "IFNG": 0,
    }

    for column, treatment, _, _ in sample_columns:
        treatment_counts[treatment] += 1
        index: int = treatment_counts[treatment]

        rename_mapping[column] = (
            f"{TREATMENT_NAMES[treatment]}{index}"
        )

    result: pl.DataFrame = (
        dataframe
        .select(selected_columns)
        .rename(rename_mapping)
    )

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    result.write_csv(output_file)

    # Preserve experiment/replicate identities for a later statistical design.
    pl.DataFrame([
        {"sample_id": rename_mapping[column], "source_sample": column,
         "condition": TREATMENT_NAMES[treatment],
         "experiment": experiment, "replicate": replicate}
        for column, treatment, experiment, replicate in sample_columns
    ]).write_csv(output_file.with_name("gse281833_sample_metadata.csv"))

    return result


if __name__ == "__main__":
    count_matrix: pl.DataFrame = prepare_count_matrix(
        source_file=SOURCE_FILE,
        output_file=OUTPUT_FILE,
    )

    print(f"Saved {count_matrix.height:,} genes to: {OUTPUT_FILE}")
