from pathlib import Path

import polars as pl

from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.manuscript.panels.fig_02.aire_distribution.confg import NEIGHBOR_DISTANCE_KB, TOP_UP_GENES, TOP_DOWN_GENES, \
    PERMUTATIONS, PERMUTATIONS_FILE, DATA_DIR
from src.manuscript.panels.fig_02.aire_distribution.shared.count_genes_with_neighbor import count_genes_with_neighbor


def main() -> None:
    # ------------------------------------------------------------------------------------------
    # PREPARE DATA
    # ------------------------------------------------------------------------------------------

    out_file: Path = PERMUTATIONS_FILE

    datasets: list[DatasetMetadata] = get_mtec_high()

    # ------------------------------------------------------------------------------------------
    # SETTINGS
    # ------------------------------------------------------------------------------------------

    permutations: int = PERMUTATIONS
    neighbor_distance_bp: int = NEIGHBOR_DISTANCE_KB * 1_000

    effect_sizes: dict[str, int] = {
        "activated": TOP_UP_GENES,
        "repressed": TOP_DOWN_GENES,
    }

    # ------------------------------------------------------------------------------------------
    # CALCULATE PERMUTATIONS
    # ------------------------------------------------------------------------------------------

    permutation_results: list[dict[str, str | int]] = []

    for dataset in datasets:
        tec_df: pl.DataFrame = dataset.get_standardized_df()

        genomic_df: pl.DataFrame = (
            tec_df
            .filter(
                pl.col("chromosome").is_not_null(),
                pl.col("start").is_not_null(),
                pl.col("stop").is_not_null(),
            )
            .select(
                "gene_symbol",
                "chromosome",
                "start",
                "stop",
            )
        )

        for aire_effect, n_genes in effect_sizes.items():
            if genomic_df.height < n_genes:
                raise ValueError(f"Dataset {dataset.dataset_id!r} has only {genomic_df.height} positioned genes, but {n_genes} are required for {aire_effect!r} permutations.")

            for permutation in range(permutations):
                random_genes_df: pl.DataFrame = genomic_df.sample(
                    n=n_genes,
                    with_replacement=False,
                    shuffle=True,
                    seed=permutation,
                )

                n_with_neighbor: int = count_genes_with_neighbor(random_genes_df, max_distance_bp=neighbor_distance_bp)

                permutation_results.append(
                    {
                        "dataset_id": dataset.dataset_id,
                        "aire_effect": aire_effect,
                        "permutation": permutation + 1,
                        "n_genes": n_genes,
                        "n_with_neighbor": n_with_neighbor,
                    }
                )

    result_df: pl.DataFrame = pl.DataFrame(permutation_results)

    # ------------------------------------------------------------------------------------------
    # SAVE
    # ------------------------------------------------------------------------------------------

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    result_df.write_csv(out_file)

if __name__ == "__main__":
    main()
