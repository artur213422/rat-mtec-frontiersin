from pathlib import Path

from config.paths import TEMP_FOLDER
from src.common.feratures.aPEAR.run_aPEAR import run_apear, APearInput
import polars as pl

from src.common.paths import PROJECT_ROOT

GO_FILE_STEM_TO_ABBREVIATION: dict[str, str] = {
    "go_biological_process_gsea": "BP",
    "go_cellular_component_gsea": "CC",
    "go_molecular_function_gsea": "MF",
}

DIRECTION_TO_GSEA_MODE: dict[str, str] = {
    "upregulated": "pos",
    "downregulated": "neg",
}


def get_standard_gsea_files(go_pathway_folder: Path) -> list[Path]:
    """Return the complete set of standard (scoreType='std') GO GSEA files."""
    expected_files: set[Path] = {
        go_pathway_folder / f"{file_stem}.csv"
        for file_stem in GO_FILE_STEM_TO_ABBREVIATION
    }
    missing_files: list[str] = sorted(
        path.name for path in expected_files if not path.is_file()
    )

    if missing_files:
        raise FileNotFoundError(
            "Missing standard GSEA inputs: "
            f"{missing_files}. Run 02_go_gsea.py first."
        )

    return sorted(expected_files)


def write_standard_gsea_direction(
    *,
    gsea_files: list[Path],
    direction: str,
    output_folder: Path,
) -> list[Path]:
    """Split standard GSEA by NES sign without recalculating NES or FDR."""
    try:
        gsea_mode: str = DIRECTION_TO_GSEA_MODE[direction]
    except KeyError as error:
        raise ValueError(f"Unsupported GSEA direction: {direction!r}") from error

    output_folder.mkdir(parents=True, exist_ok=True)
    directional_files: list[Path] = []

    for gsea_file in gsea_files:
        gsea_df: pl.DataFrame = pl.read_csv(gsea_file)
        required_columns: set[str] = {"dataset_id", "term", "nes", "fdr_q_value"}
        missing_columns: set[str] = required_columns - set(gsea_df.columns)
        if missing_columns:
            raise ValueError(
                f"Standard GSEA input {gsea_file} is missing columns: "
                f"{sorted(missing_columns)}"
            )

        nes_filter = pl.col("nes") > 0 if direction == "upregulated" else pl.col("nes") < 0
        directional_df: pl.DataFrame = gsea_df.filter(
            pl.col("nes").is_not_null()
            & pl.col("nes").is_finite()
            & nes_filter
        )

        output_file: Path = output_folder / f"{gsea_file.stem}_{gsea_mode}.csv"
        directional_df.write_csv(output_file)
        directional_files.append(output_file)

    return directional_files

def get_go_abbreviation(gsea_file: Path, *, gsea_mode: str) -> str:
    file_stem: str = gsea_file.stem.removesuffix(f"_{gsea_mode}")
    try:
        return GO_FILE_STEM_TO_ABBREVIATION[file_stem]
    except KeyError as error:
        raise ValueError(f"Unrecognised GO GSEA input file: {gsea_file}") from error


def validate_gsea_files(gsea_files: list[Path], *, gsea_mode: str, go_pathway_folder: Path) -> None:
    expected_files: set[Path] = {
        go_pathway_folder / f"{file_stem}_{gsea_mode}.csv"
        for file_stem in GO_FILE_STEM_TO_ABBREVIATION
    }
    actual_files: set[Path] = set(gsea_files)

    if actual_files != expected_files:
        missing_files: list[str] = sorted(path.name for path in expected_files - actual_files)
        unexpected_files: list[str] = sorted(path.name for path in actual_files - expected_files)
        raise FileNotFoundError(f"Invalid GO GSEA inputs for mode={gsea_mode!r}. Missing: {missing_files}; unexpected: {unexpected_files}. Run 03_run_one_sided_gsea.py first.")


def run_direction(*, direction: str, gsea_mode: str, gsea_files: list[Path], apear_input_folder: Path) -> None:
    output_folder: Path = apear_input_folder / direction
    output_folder.mkdir(parents=True, exist_ok=True)

    for gsea_file in gsea_files:
        go_abbreviation: str = get_go_abbreviation(gsea_file, gsea_mode=gsea_mode,)
        gsea_df: pl.DataFrame = pl.read_csv(gsea_file)

        if "dataset_id" not in gsea_df.columns:
            raise ValueError(f"GSEA input is missing the dataset_id column: {gsea_file}")
        if gsea_df.get_column("dataset_id").null_count() > 0:
            raise ValueError(f"GSEA input contains a null dataset_id: {gsea_file}")

        for (dataset_id_value,), dataset_df in gsea_df.group_by("dataset_id", maintain_order=True,):
            dataset_id: str = str(dataset_id_value)
            output_prefix: str = f"{dataset_id}_{go_abbreviation}"
            dataset_gsea_file: Path = PROJECT_ROOT / TEMP_FOLDER / f"{output_prefix}_{direction}.csv"


            # aPEAR requires one dataset per input file.
            dataset_df.write_csv(dataset_gsea_file)

            run_apear(
                APearInput(
                    gsea_file=dataset_gsea_file,
                    output_svg=output_folder / f"{output_prefix}.svg",
                    output_clusters=(
                        output_folder / f"{output_prefix}_clusters.csv"
                    ),
                )
            )
