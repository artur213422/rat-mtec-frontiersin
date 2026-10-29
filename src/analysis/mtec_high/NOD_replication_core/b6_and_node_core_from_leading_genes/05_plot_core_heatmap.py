"""Plot leading-edge overlap classes; BALB/c remains in the annotated CSV only."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from src.common.paths import PROJECT_ROOT

output_folder: Path = PROJECT_ROOT / "data/mtec_high/derived_data/NOD_replication_core/b6_and_nod_core_from_leading_genes"
df = pl.read_csv(output_folder / "standardized/nod_and_b6_core.csv", schema_overrides={"ortholog_id": pl.String})
datasets = [
    ("mouse_gse33878_mtec_high", "NOD\nGSE33878"),
    ("mouse_gse8564_nod_mtec_high", "NOD\nGSE8564"),
    ("mouse_gse85_mtec_high", "Other\nGSE85"),
    ("mouse_gse8564_b6_mtec_high", "B6\nGSE8564"),
    ("mouse_gse2585_mtec_high", "Other\nGSE2585"),
    ("rat_mtec_high", "Rat"),
]
classes = [("nod_only", "NOD-only genes"), ("shared", "Shared genes"), ("b6_rat_only", "B6/Rat-only genes")]
lookup = {(r["leading_gene"], r["dataset_id"]): r["log2_fold_change"] for r in df.iter_rows(named=True)}
panels = []
order_rows = []
for group, title in classes:
    genes = sorted(df.filter(pl.col("overlap_class") == group)["leading_gene"].unique().to_list())
    # Sort by mean NOD effect, then symbol, so row order is reproducible.
    def sort_key(gene):
        values = [lookup.get((gene, dataset)) for dataset, _ in datasets[:2]]
        values = [v for v in values if v is not None and np.isfinite(v)]
        return (float(np.mean(values)) if values else float("inf"), gene)
    genes.sort(key=sort_key)
    matrix = np.array([[lookup.get((gene, dataset)) for dataset, _ in datasets] for gene in genes], dtype=float).reshape(-1, len(datasets))
    panels.append((title, genes, matrix))
    order_rows.extend({"overlap_class": group, "row": i + 1, "leading_gene": gene} for i, gene in enumerate(genes))

fig, axes = plt.subplots(3, 1, figsize=(12, 19), gridspec_kw={"height_ratios": [max(len(g), 10) for _, g, _ in panels]})
cmap = plt.get_cmap("RdBu_r").copy()
cmap.set_bad("#d9d9d9")
for ax, (title, genes, matrix) in zip(axes, panels):
    displayed = matrix if len(genes) else np.full((1, len(datasets)), np.nan)
    im = ax.imshow(np.ma.masked_invalid(displayed), aspect="auto", interpolation="nearest", cmap=cmap, vmin=-0.6, vmax=0.6)
    ax.set_title(f"{title} (n={len(genes)})", loc="left", weight="bold", fontsize=15)
    ax.axvline(1.5, color="black", linewidth=2)
    ax.set_yticks([])
    ax.set_xticks([])
axes[-1].set_xticks(range(len(datasets)), [label for _, label in datasets], fontsize=11)
axes[0].text(1/6, 1.09, "NOD datasets", transform=axes[0].transAxes, ha="center", fontsize=13)
axes[0].text(2/3, 1.09, "Other / B6 / Rat datasets", transform=axes[0].transAxes, ha="center", fontsize=13)
fig.suptitle("Core pathway genes grouped by overlap class\nCell values: log2 fold change", fontsize=18, weight="bold", y=0.985)
fig.subplots_adjust(left=0.04, right=0.86, top=0.91, bottom=0.06, hspace=0.24)
cax = fig.add_axes([0.89, 0.40, 0.02, 0.20])
fig.colorbar(im, cax=cax, label="log2 fold change", extend="both")
fig.text(0.04, 0.015, "Grey: unavailable measurement. Colours saturate at ±0.6. Membership does not imply gene-level significance.", fontsize=9)
figures = output_folder / "figures"
figures.mkdir(parents=True, exist_ok=True)
fig.savefig(figures / "04_core_heatmap.png", dpi=180)
fig.savefig(figures / "04_core_heatmap.pdf")
plt.close(fig)
pl.DataFrame(order_rows, schema={"overlap_class": pl.String, "row": pl.Int64, "leading_gene": pl.String}).write_csv(figures / "04_heatmap_row_order.csv")
print(f"Saved heatmap and row order to {figures}")
