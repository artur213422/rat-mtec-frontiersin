"""Plot median response groups using pathway-table labels; BALB/c is comparison-only."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from src.common.paths import PROJECT_ROOT

output_folder: Path = PROJECT_ROOT / "data/mtec_high/derived_data/NOD_replication_core/b6_and_nod_core_from_leading_genes"
df = pl.read_csv(output_folder / "standardized/grouped_core_genes.csv", schema_overrides={"ortholog_id": pl.String})
datasets = [
    ("mouse_gse33878_mtec_high", "NOD I\nGSE33878"),
    ("mouse_gse8564_nod_mtec_high", "NOD II\nGSE8564"),
    ("mouse_gse85_mtec_high", "B6 I\nGSE85"),
    ("mouse_gse8564_b6_mtec_high", "B6 II\nGSE8564"),
    ("mouse_gse2585_mtec_high", "B6 III\nGSE2585"),
    ("rat_mtec_high", "Rat"),
    ("mouse_gse8564_balbc_mtec_high", "BALB/c\nGSE8564"),
]
classes = [("both_down", "Both group medians down"),
           ("opposite_directions", "Opposite group median directions"),
           ("both_up", "Both group medians up")]
if df.filter(pl.col("response_group") == "unclassified").height:
    classes.append(("unclassified", "Unclassified: zero or missing median"))
if df.select("leading_gene", "dataset_id").is_duplicated().any():
    raise ValueError("Expected one measurement per gene and dataset.")
lookup = {(r["leading_gene"], r["dataset_id"]): r["log2_fold_change"] for r in df.iter_rows(named=True)}
panels = []
order_rows = []
for group, title in classes:
    genes = sorted(df.filter(pl.col("response_group") == group)["leading_gene"].unique().to_list())
    # Sort by table label NOD I (actually B6 x NOD F1), then symbol; missing last.
    def sort_key(gene):
        value = lookup.get((gene, datasets[0][0]))
        return (float(value) if value is not None and np.isfinite(value) else float("inf"), gene)
    genes.sort(key=sort_key)
    matrix = np.array([[lookup.get((gene, dataset)) for dataset, _ in datasets] for gene in genes], dtype=float).reshape(-1, len(datasets))
    panels.append((title, genes, matrix))
    order_rows.extend({"response_group": group, "row": i + 1, "leading_gene": gene} for i, gene in enumerate(genes))

fig, axes = plt.subplots(len(panels), 1, figsize=(12, 19), gridspec_kw={"height_ratios": [max(len(g), 10) for _, g, _ in panels]})
cmap = plt.get_cmap("RdBu_r").copy()
cmap.set_bad("#d9d9d9")
for ax, (title, genes, matrix) in zip(axes, panels):
    displayed = matrix if len(genes) else np.full((1, len(datasets)), np.nan)
    im = ax.imshow(np.ma.masked_invalid(displayed), aspect="auto", interpolation="nearest", cmap=cmap, vmin=-0.6, vmax=0.6)
    ax.set_title(f"{title} (n={len(genes)})", loc="left", weight="bold", fontsize=15)
    ax.axvline(1.5, color="black", linewidth=2)
    ax.axvline(5.5, color="black", linewidth=2)
    ax.set_yticks([])
    ax.set_xticks([])
axes[-1].set_xticks(range(len(datasets)), [label for _, label in datasets], fontsize=11)
axes[0].text(1/7, 1.09, "NOD group", transform=axes[0].transAxes, ha="center", fontsize=13)
axes[0].text(4/7, 1.09, "B6 mouse / Rat group", transform=axes[0].transAxes, ha="center", fontsize=13)
axes[0].text(6.5/7, 1.09, "Other", transform=axes[0].transAxes, ha="center", fontsize=11)
fig.suptitle("Core pathway genes grouped by median response\nCell values: log2 fold change", fontsize=18, weight="bold", y=0.985)
fig.subplots_adjust(left=0.04, right=0.86, top=0.91, bottom=0.08, hspace=0.24)
cax = fig.add_axes([0.89, 0.40, 0.02, 0.20])
fig.colorbar(im, cax=cax, label="log2 fold change", extend="both")
fig.text(0.04, 0.028, "Grey: missing. Scale: ±0.6. Rows ordered by NOD I. Groups use median signs, not significance; BALB/c excluded.", fontsize=9)
fig.text(0.04, 0.015, "Table labels retained: NOD I = B6×NOD F1; B6 I = B6/129 F2; B6 III = B6 WT vs mixed-background KO.", fontsize=9)
figures = output_folder / "figures"
figures.mkdir(parents=True, exist_ok=True)
fig.savefig(figures / "06_grouped_heatmap.png", dpi=180)
fig.savefig(figures / "06_grouped_heatmap.pdf")
plt.close(fig)
pl.DataFrame(order_rows, schema={"response_group": pl.String, "row": pl.Int64, "leading_gene": pl.String}).write_csv(figures / "06_grouped_heatmap_row_order.csv")
print(f"Saved heatmap and row order to {figures}")
