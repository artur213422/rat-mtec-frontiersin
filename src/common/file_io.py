import gzip
import json
import shutil
from typing import Any, TextIO
from pathlib import Path
import pandas as pd
import polars as pl
import requests

GAF_COLUMNS = [
    "db",
    "db_object_id",
    "gene_symbol",
    "qualifier",
    "go_term",
    "reference",
    "evidence_code",
    "with_from",
    "aspect",
    "gene_name",
    "gene_synonym",
    "object_type",
    "taxon",
    "date",
    "assigned_by",
    "annotation_extension",
    "gene_product_form_id",
]

def read_gaf_gz_file(file: Path) -> pl.DataFrame:
    return pl.read_csv(
        file,
        separator="\t",
        has_header=False,
        comment_prefix="!",
        new_columns=GAF_COLUMNS,
        schema_overrides={column: pl.Utf8 for column in GAF_COLUMNS},
    )


def download_from_web(url: str, output_path: Path):
    headers = {
        "User-Agent": "Mozilla/5.0"
    }

    with requests.get(url, headers=headers, stream=True, timeout=60) as r:
        r.raise_for_status()

        with output_path.open("wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)

def unpack_file(file_path: Path, output_dir: Path | None = None) -> Path:
    if not file_path.exists():
        raise FileNotFoundError(file_path)

    if output_dir is None:
        output_dir = file_path.parent / file_path.stem

    output_dir.mkdir(parents=True, exist_ok=True)

    if file_path.suffix == ".gz" and not file_path.name.endswith(".tar.gz"):
        output_path: Path = output_dir / file_path.with_suffix("").name

        with gzip.open(file_path, "rb") as source:
            with output_path.open("wb") as target:
                shutil.copyfileobj(source, target)

        print(f"Extracted {file_path.name} -> {output_path}")
        return output_path

    shutil.unpack_archive(
        filename=file_path,
        extract_dir=output_dir,
    )

    print(f"Extracted {file_path.name} -> {output_dir}")
    return output_dir

def delete_file(file_path: Path):
    if not file_path.exists():
        return

    if not file_path.is_file():
        raise ValueError(f"Expected a file, got: {file_path}")

    file_path.unlink()

    print(f"Deleted: {file_path}")

def load_json(file_path: Path) ->  dict[str, Any]:
    try:
        with file_path.open("r", encoding="utf-8") as file:
            data: Any = json.load(file)


    except FileNotFoundError as error:

        raise FileNotFoundError(f"JSON file was not found: {file_path}") from error


    except json.JSONDecodeError as error:

        raise ValueError(

            f"File is not valid JSON: {file_path}. "

            f"Error at line {error.lineno}, column {error.colno}: {error.msg}"

        ) from error

    if not isinstance(data, dict):
        raise TypeError(
            f"Expected JSON top-level object to be dict[str, Any], "
            f"but got {type(data).__name__}: {file_path}"
        )

    return data


def save_csv_to_excel(file_path: Path, output_file: Path | None = None) -> Path:

    if output_file is None:
        output_file = file_path.with_suffix(".xlsx")

    output_file.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(file_path)

    df.to_excel(
        output_file,
        index=False,
        engine="openpyxl",
    )

    print(f"Saved Excel file: {output_file}")
    return output_file

def _open_text(path: Path) -> TextIO:
    """Open plain text or .gz GEO files."""
    if path.suffix == ".gz":
        return gzip.open(path, mode="rt", encoding="utf-8", errors="replace")
    return path.open(mode="rt", encoding="utf-8", errors="replace")




def read_gmt_file(file: Path) -> pl.DataFrame:
    rows: list[dict[str, str]] = []

    with file.open("r", encoding="utf-8") as handle:
        for line in handle:
            columns: list[str] = line.rstrip().split("\t")

            if len(columns) < 3:
                continue

            pathway: str = columns[0]
            genes: list[str] = columns[2:]

            for gene_symbol in genes:
                if not gene_symbol:
                    continue

                rows.append(
                    {
                        "pathway": pathway,
                        "gene_symbol": gene_symbol,
                    }
                )

    return pl.DataFrame(
        rows,
        schema={
            "pathway": pl.String,
            "gene_symbol": pl.String,
        },
    )

def ask_yes_no(prompt: str) -> bool:
    while True:
        answer: str = input(f"{prompt} [y/n]: ").strip().lower()

        if answer in {"y", "yes"}:
            return True

        if answer in {"n", "no"}:
            return False

        print("Please answer 'y' or 'n'.")


def delete_all_files(folder_path: Path, ask_permission: bool = False) -> None:
    if not folder_path.exists():
        raise FileNotFoundError(f"Folder does not exist: {folder_path}")

    if not folder_path.is_dir():
        raise NotADirectoryError(f"Path is not a directory: {folder_path}")

    if ask_permission and not ask_yes_no(
            f"Are you sure you want to delete all files in this directory?\n"
            f"{folder_path}"
    ):
        return

    for path in folder_path.iterdir():
        if path.is_file():
            path.unlink()