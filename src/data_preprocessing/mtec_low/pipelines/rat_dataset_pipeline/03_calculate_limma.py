from pathlib import Path

from src.common.dataset_catalog.api.get_datasets import get_rat_tec_low
from src.common.dataset_catalog.internals.types import DatasetMetadata

import polars as pl

from src.common.feratures.limma.run_limma import run_limma, LimmaInput
from src.common.paths import resolve_dir_from_root

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

rat_dataset: DatasetMetadata = get_rat_tec_low()

tec_df: pl.DataFrame = rat_dataset.get_standardized_df()

output_file: Path = resolve_dir_from_root("data/mtec_low/standardized/rat") / f"{rat_dataset.dataset_id}.csv"
##########################################################################################
#                                           CALCULATE
##########################################################################################

limma_input: LimmaInput = LimmaInput(
    frame=tec_df,
    output_file=output_file,
    data_type="rma",
    merge=True
)

tec_df = run_limma(limma_input)


##########################################################################################
#                                           SAVE
##########################################################################################


tec_df.write_csv(output_file)
print("Done")