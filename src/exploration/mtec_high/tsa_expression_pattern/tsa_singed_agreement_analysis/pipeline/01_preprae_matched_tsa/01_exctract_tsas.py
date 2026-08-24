from pathlib import Path
from statistics import mean

import polars as pl

from src.common.dataframe_ops import get_genotype_expression_values
from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.paths import try_find_file_in_project, resolve_dir_from_root
from src.common.types.enums import Genotypes

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

ranked_tsa_file: Path = try_find_file_in_project("ranked_tsa.csv")
ranked_tsa_df: pl.DataFrame = pl.read_csv(ranked_tsa_file)
tsa_file: Path = resolve_dir_from_root("data/mtec_high/analysis/tsa_expression_pattern/tsa_singed_agreement_analysis") / "defined_tsa.csv"


##########################################################################################
#                                       SETTINGS
##########################################################################################

group_max_log2: int | None = 6         # the gene is actually expressed reasonably strongly somewhere in the tissue/group.
mean_enrichment_log2: int = 1   # the group is on average at least 2-fold higher than the rest of the body.
max_specificity_log2: int = 1   # even the strongest outside tissue is at least 2-fold lower than the strongest tissue inside the group.


##########################################################################################
#                                       PROCESS DATA
##########################################################################################
filtered_ranked_tsa_df: pl.DataFrame = ranked_tsa_df.filter(
    (pl.col("mean_enrichment_log2") >= mean_enrichment_log2)
     & (pl.col("max_specificity_log2") >= max_specificity_log2)
)

if group_max_log2 is not None:
    filtered_ranked_tsa_df: pl.DataFrame = ranked_tsa_df.filter(
        pl.col("group_max_log2") >= 1
    )

##########################################################################################
#                                           SAVE
##########################################################################################

filtered_ranked_tsa_df.write_csv(tsa_file)
