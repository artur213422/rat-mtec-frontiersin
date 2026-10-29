from zipfile import Path

import polars as pl

from src.common.paths import find_file_in_project

chromosome_file: Path = find_file_in_project("rat_chromosome_data.txt")

def attach_chromosome_data(tec_df: pl.DataFrame) -> pl.DataFrame:
    """Attach rat transcript-cluster genomic coordinates by array ``ID``."""
    columns_to_attach = ["chromosome", "strand", "start", "stop"]

    if "ID" not in tec_df.columns:
        raise ValueError("tec_df must contain an 'ID' column.")

    chromosome_df = pl.read_csv(
        chromosome_file,
        separator="\t",
        columns=["ID", *columns_to_attach],
    )

    duplicate_ids = chromosome_df.group_by("ID").len().filter(pl.col("len") > 1)
    if duplicate_ids.height:
        raise ValueError("rat_chromosome_data.txt contains duplicate ID values.")

    return tec_df.join(chromosome_df, on="ID", how="left")