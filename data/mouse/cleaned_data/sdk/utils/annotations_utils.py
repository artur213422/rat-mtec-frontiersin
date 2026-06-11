from io import StringIO
from pathlib import Path

import polars as pl
from polars import DataFrame

def get_gene_symbol(probe_id: str, annot: DataFrame, selected_col: str) -> str | None:
    result = (
        annot
        .filter(pl.col("ID") == probe_id)
        .select(selected_col)
    )


    if result.height == 0:
        return None

    gene_symbol: str | None = result.item()

    if gene_symbol is None or gene_symbol == "":
        return None

    return gene_symbol


def get_platform_id(file: str) -> str:
    platform_id: str | None = None

    with open(file, encoding="utf-8") as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith("!Sample_platform_id"):
            parts = line.rstrip("\n").split("\t")
            platform_id = parts[1]
            break

    if platform_id is None:
        raise ValueError(f"Failed to extract platform from file={file}")

    del lines

    return platform_id.strip('"')

def get_matrix(file: str) -> DataFrame:
    with open(file, encoding="utf-8") as f:
        lines = f.readlines()

    sample_titles: list[str] | None = None

    for line in lines:
        if line.startswith("!Sample_title"):
            parts = line.rstrip("\n").split("\t")
            sample_titles = [x.strip('"') for x in parts[1:]]
            break


    if sample_titles is None:
        raise ValueError(f"Failed to extract sample titles from file={file}")

    table_start = None

    for i, line in enumerate(lines):
        if line.startswith("!series_matrix_table_begin"):
            table_start = i + 1
            break

    if table_start is None:
        raise ValueError(f"Failed to find matrix section in {file}")

    table_end: int | None = None

    for i, line in enumerate(lines):
        if line.startswith("!series_matrix_table_end"):
            table_end = i
            break

    if table_end is None:
        raise ValueError(f"Failed to find matrix end in {file}")

    matrix_text = "".join(lines[table_start:table_end])

    # Free memory
    del lines

    # Read matrix
    df = pl.read_csv(
        StringIO(matrix_text),
        separator="\t",
        quote_char='"'
    )

    if df.is_empty():
        raise ValueError(f"No expression data found in {file}")

    # Rename GSM columns
    gsm_cols = df.columns[1:]

    rename_map = {
        gsm: title
        for gsm, title in zip(gsm_cols, sample_titles)
    }

    df = df.rename(rename_map)

    return df

def get_annotation_table(folder: str, platform_id: str) -> DataFrame:
    folder_path = Path(folder)

    matches = list(folder_path.glob(f"{platform_id}-*"))

    if len(matches) == 0:
        raise FileNotFoundError(f"No annotation file found for platform_id={platform_id} in folder={folder}")

    if len(matches) > 1:
        raise ValueError(
            f"Multiple annotation files found for platform_id={platform_id}: "
            f"{[f.name for f in matches]}"
        )

    annot_file = matches[0]

    annot = pl.read_csv(
        annot_file,
        separator="\t",
        comment_prefix="#",
        infer_schema_length=0,
    )

    if annot.is_empty():
        raise ValueError(f"Annotation table is empty: {annot_file}")

    return annot


def fix_GSE33878(df: DataFrame) -> DataFrame:
    cleaned_symbols: list[str | None] = []

    for gene_symbol in df["gene_symbol"]:
        if gene_symbol is None:
            cleaned_symbols.append(None)
            continue

        parts = [part.strip() for part in gene_symbol.split("//")]

        if len(parts) < 2:
            cleaned_symbols.append(None)
            continue

        cleaned_symbols.append(parts[1])

    return df.with_columns(
        pl.Series("gene_symbol", cleaned_symbols)
    )