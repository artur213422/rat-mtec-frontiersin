from pathlib import Path

import polars as pl

from src.manuscript.panels.fig_02.aire_distribution.confg import NEIGHBOR_COUNTS_FILE, TOP_GENES_FILE, \
    NEIGHBOR_DISTANCE_KB, DATA_DIR
from src.manuscript.panels.fig_02.aire_distribution.shared.count_genes_with_neighbor import count_genes_with_neighbor


def main() -> None:
    # ------------------------------------------------------------------------------------------
    # PREPARE DATA
    # ------------------------------------------------------------------------------------------

    top_genes: Path = TOP_GENES_FILE

    out_file: Path = NEIGHBOR_COUNTS_FILE

    # Keep counts and denominators restricted to genes with usable coordinates.
    top_df: pl.DataFrame = pl.read_csv(top_genes).drop_nulls(
        subset=["chromosome", "start", "stop"],
    )

    neighbor_distance_bp: int = NEIGHBOR_DISTANCE_KB * 1_000

    # ------------------------------------------------------------------------------------------
    # COUNT NEIGHBORS
    # ------------------------------------------------------------------------------------------

    counts: list[dict[str, str | int | float]] = []

    for (dataset_id, aire_effect), group_df in top_df.group_by(["dataset_id", "aire_effect"], maintain_order=True,):
        n_genes: int = group_df.height
        n_with_neighbor: int = count_genes_with_neighbor(
            group_df, max_distance_bp=neighbor_distance_bp,
        ) if n_genes >= 2 else 0

        counts.append(
            {
                "dataset_id": dataset_id,
                "aire_effect": aire_effect,
                "n_genes": n_genes,
                "n_with_neighbor": n_with_neighbor,
                "fraction_with_neighbor": n_with_neighbor / n_genes,
            }
        )

    result_df: pl.DataFrame = pl.DataFrame(counts)

    # ------------------------------------------------------------------------------------------
    # SAVE
    # ------------------------------------------------------------------------------------------

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    result_df.write_csv(out_file)

if __name__ == "__main__":
    main()
