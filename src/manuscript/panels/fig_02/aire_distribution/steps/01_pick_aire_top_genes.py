from pathlib import Path
import polars as pl

from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.manuscript.panels.fig_02.aire_distribution.confg import TOP_GENES_FILE, TOP_DOWN_GENES, TOP_UP_GENES, \
    LOG2_FC_THRESHOLD, DATA_DIR


def main() -> None:
    # ----------------------------------------------------------------------------------------------
    ## PREPARE DATA
    # ----------------------------------------------------------------------------------------------

    out_file: Path = TOP_GENES_FILE

    datasets: list[DatasetMetadata] = get_mtec_high()

    # ----------------------------------------------------------------------------------------------
    ## SETTINGS
    # ----------------------------------------------------------------------------------------------

    top_up_genes: int = TOP_UP_GENES
    top_down_genes: int = TOP_DOWN_GENES

    p_val_threshold: float = 0.05
    log2_fc_threshold: float = LOG2_FC_THRESHOLD

    selected_columns: list[str] = [
        "gene_symbol",
        "chromosome",
        "strand",
        "start",
        "stop",
        "log2_fold_change",
        "t_statistic",
        "p_val",
        "p_adj",
        "dataset_id",
        "aire_effect",
    ]

    # ----------------------------------------------------------------------------------------------
    ## PICK TOP
    # ----------------------------------------------------------------------------------------------

    picked: list[pl.DataFrame] = []

    for dataset in datasets:

        # Apply the same coordinate eligibility as the random background before ranking.
        tec_df: pl.DataFrame = dataset.get_standardized_df().drop_nulls(
            subset=["chromosome", "start", "stop"],
        )

        significant_df: pl.DataFrame = tec_df.filter(
            (pl.col("p_val") < p_val_threshold)
            & (pl.col("log2_fold_change").abs() >= log2_fc_threshold)
        )

        aire_activated_df: pl.DataFrame = significant_df.filter(pl.col("log2_fold_change") < 0)

        aire_repressed_df: pl.DataFrame = significant_df.filter(pl.col("log2_fold_change") > 0)

        if aire_activated_df.height < top_up_genes:
            raise ValueError(f"Dataset {dataset.dataset_id!r} has only {aire_activated_df.height} positioned Aire-activated DEGs, but {top_up_genes} are required.")

        if aire_repressed_df.height < top_down_genes:
            raise ValueError(f"Dataset {dataset.dataset_id!r} has only {aire_repressed_df.height} positioned Aire-repressed DEGs, but {top_down_genes} are required.")

        aire_activated: pl.DataFrame = (
            aire_activated_df
            .sort("log2_fold_change")
            .head(top_up_genes)
            .with_columns(
                pl.lit(dataset.dataset_id).alias("dataset_id"),
                pl.lit("activated").alias("aire_effect"),
            )
            .select(selected_columns)
        )

        aire_repressed: pl.DataFrame = (
            aire_repressed_df
            .sort("log2_fold_change", descending=True)
            .head(top_down_genes)
            .with_columns(
                pl.lit(dataset.dataset_id).alias("dataset_id"),
                pl.lit("repressed").alias("aire_effect"),
            )
            .select(selected_columns)
        )

        picked.append(pl.concat([aire_activated, aire_repressed], how="vertical"))

    result_df: pl.DataFrame = pl.concat(picked, how="vertical",)

    # ----------------------------------------------------------------------------------------------
    ## SAVE
    # ----------------------------------------------------------------------------------------------

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    result_df.write_csv(out_file)

if __name__ == "__main__":
    main()
