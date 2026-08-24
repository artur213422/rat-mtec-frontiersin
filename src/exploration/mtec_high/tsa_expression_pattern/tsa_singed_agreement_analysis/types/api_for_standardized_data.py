from enum import StrEnum
from pathlib import Path

import polars as pl

from src.common.paths import find_file_in_project

class AireStatus(StrEnum):
    DOWNREGULATED = "downregulated"
    NOT_CHANGED = "not_changed"
    UPREGULATED = "upregulated"


class AireGlobalStatus(StrEnum):
    CHANGED = "changed"
    UNCHANGED = "unchanged"


AIRE_STATUS_TYPE: pl.Enum = pl.Enum(AireStatus)
AIRE_GLOBAL_STATUS_TYPE: pl.Enum = pl.Enum(AireGlobalStatus)

def get_all_tsa_df() -> pl.DataFrame:
    file: Path = find_file_in_project("matched_tsa_list.csv")

    return pl.read_csv(file)

def get_tsa_changed_in_any_dataset_df() -> pl.DataFrame:
    tsa_df: pl.DataFrame = get_all_tsa_df()

    return tsa_df.filter(
        pl.col("aire_status_global") == AireGlobalStatus.CHANGED
    )