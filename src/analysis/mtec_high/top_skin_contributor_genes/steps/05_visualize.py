"""Rank contributors by fold change and compare all datasets using clean orthologs."""
from pathlib import Path
import math
import re

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import polars as pl

from src.analysis.mtec_high.top_skin_contributor_genes.config import (
    CONTRIBUTOR_P_VAL_THRESHOLD,
)
from src.common.dataset_catalog.api.get_datasets import (
    get_mouse_tec_high, get_rat_tec_high, read_file_with_ortholog_id,
)
from src.common.feratures.visualization.heatmap.heatmap import Heatmap
from src.common.feratures.visualization.heatmap.types.heatmap_data import HeatmapData
from src.common.paths import PROJECT_ROOT

BASE = PROJECT_ROOT / "data/mtec_high/derived_data/top_skin_contributor_genes"
TOP_GENES = 30


def select_top(candidates: pl.DataFrame) -> pl.DataFrame:
    """Mean absolute fold change over supporting datasets; no t-statistic ranking."""
    return (
        candidates.filter(pl.col("log2_fold_change").is_finite())
        .unique(subset=["dataset_id", "gene_symbol"])
        .group_by("gene_symbol")
        .agg(
            pl.col("log2_fold_change").abs().mean().alias("ranking_score"),
            pl.col("dataset_id").n_unique().alias("supporting_datasets"),
            pl.col("ortholog_id").drop_nulls().first(),
        )
        .sort(["ranking_score", "gene_symbol"], descending=[True, False])
        .head(TOP_GENES).with_row_index("rank", offset=1)
    )


def comparison_data(selected: pl.DataFrame, reference: pl.DataFrame,
                    dataset_ids: list[str], focal: pl.DataFrame) -> pl.DataFrame:
    selected = selected.with_columns(
        (pl.col("rank").cast(pl.String).str.pad_start(2, "0")
         + pl.lit(". ") + pl.col("gene_symbol")).alias("row_label")
    )
    grid = selected.join(pl.DataFrame({"dataset_id": dataset_ids}), how="cross")
    grid = grid.join(
        reference.filter(pl.col("has_clean_ortholog"))
        .select("dataset_id", "ortholog_id", pl.col("gene_symbol").alias("matched_gene_symbol"),
                "log2_fold_change", "p_val"),
        on=["dataset_id", "ortholog_id"], how="left", validate="m:1",
    )
    # Within the focal species, use symbols to recover every measured dataset,
    # including non-contributors and genes without a clean cross-species mapping.
    is_rat = focal["dataset_id"][0].startswith("rat")
    same_species = reference.filter(pl.col("dataset_id").str.starts_with("rat") == is_rat)
    return grid.join(
        same_species.select(
            "dataset_id", "gene_symbol",
            pl.col("log2_fold_change").alias("focal_fc"),
            pl.col("p_val").alias("focal_p_val"),
            pl.lit(True).alias("focal_match"),
        )
        .unique(subset=["dataset_id", "gene_symbol"]),
        on=["dataset_id", "gene_symbol"], how="left", validate="m:1",
    ).with_columns(
        pl.when(pl.col("focal_match").fill_null(False)).then(pl.col("focal_fc"))
        .otherwise(pl.col("log2_fold_change")).alias("log2_fold_change"),
        pl.when(pl.col("focal_match").fill_null(False)).then(pl.col("focal_p_val"))
        .otherwise(pl.col("p_val")).alias("p_val"),
    ).drop("focal_fc", "focal_p_val", "focal_match")


def draw_comparison(grid: pl.DataFrame, dataset_ids: list[str], source: str,
                    pathway: str, folder: Path, mode: str) -> None:
    n_genes = grid["gene_symbol"].n_unique()
    fig, ax = plt.subplots(figsize=(12, max(5.5, n_genes * 0.34 + 2)))
    heatmap = Heatmap(ax, HeatmapData(grid, "log2_fold_change", "dataset_id", "row_label"))
    heatmap.show_col_labels(rotation=30, labels={
        dataset: dataset.replace("_mtec_high", "").replace("mouse_", "Mouse ").replace("rat", "Rat")
        for dataset in dataset_ids
    })
    finite = grid.filter(pl.col("log2_fold_change").is_finite())
    limit = max(1, math.ceil(finite["log2_fold_change"].abs().max()))
    heatmap.configure_color_scale(vmin=-limit, vmax=limit, tick=max(1, math.ceil(limit / 4)))
    heatmap.show_color_bar("log2 fold change", size=0.35, label_distance=35)
    heatmap.show_na_legend()
    heatmap.show_significance(
        value_threshold=1,
        p_value_col="p_val",
        p_threshold=CONTRIBUTOR_P_VAL_THRESHOLD,
    )
    heatmap.show_significance_legend()
    ax.set_title(f"{source}: {pathway}", pad=18, fontsize=14)
    fig.text(0.13, 0.035,
             "Rank: mean |log2 FC| across supporting mouse datasets" if mode == "best_mouse_core"
             else "Rank: rat |log2 FC|", fontsize=10)
    fig.subplots_adjust(left=0.18, right=0.81, bottom=0.09, top=0.76)
    folder.mkdir(parents=True, exist_ok=True)
    filename = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", f"{source}__{pathway}").rstrip(" .")
    fig.savefig(folder / f"{filename}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    grid.with_columns(pl.lit(source).alias("source"), pl.lit(pathway).alias("pathway")).write_csv(folder / f"{filename}.csv")
    print(f"Saved {mode}: {source}: {pathway} ({n_genes} genes)")


def main() -> None:
    mice = get_mouse_tec_high()
    rat = get_rat_tec_high()
    dataset_ids = [rat.dataset_id, *(m.dataset_id for m in mice)]
    reference = pl.concat([
        metadata.get_standardized_df().select(
            "gene_symbol", "ortholog_id", "has_clean_ortholog", "log2_fold_change", "p_val"
        ).with_columns(pl.lit(metadata.dataset_id).alias("dataset_id"))
        for metadata in [rat, *mice]
    ])
    mouse_candidates = read_file_with_ortholog_id(BASE / "standardized/top_skin_contributors.csv")
    mouse_candidates = mouse_candidates.filter(pl.col("dataset_id").is_in([m.dataset_id for m in mice]))
    # Step 03 removes rats. Rank all their step 02 candidates; mark significance
    # with asterisks in the heatmap instead of filtering genes before ranking.
    rat_candidates = read_file_with_ortholog_id(BASE / "source/02_annotated_leading_genes.csv").filter(
        pl.col("dataset_id") == rat.dataset_id
    )
    for mode, candidates in [("best_rat", rat_candidates), ("best_mouse_core", mouse_candidates)]:
        candidates = candidates.drop("ortholog_id", strict=False).join(
            reference.select("dataset_id", "gene_symbol", "ortholog_id", "has_clean_ortholog"),
            on=["dataset_id", "gene_symbol"], how="left", validate="m:1",
        ).with_columns(
            pl.when(pl.col("has_clean_ortholog")).then(pl.col("ortholog_id"))
            .otherwise(None).alias("ortholog_id")
        )
        for group in candidates.partition_by(["source", "pathway"], maintain_order=True):
            source, pathway = group.select("source", "pathway").row(0)
            selected = select_top(group)
            if selected.is_empty():
                print(f"Skipped {mode}: {source}: {pathway} (no finite fold changes)")
                continue
            grid = comparison_data(selected, reference, dataset_ids, group)
            draw_comparison(grid, dataset_ids, source, pathway,
                            BASE / "standardized/visualization" / mode, mode)


if __name__ == "__main__":
    main()
