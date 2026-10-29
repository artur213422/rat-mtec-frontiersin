"""Mouse genomic clusters using the same plotting style as Figure 2B."""
from pathlib import Path
import polars as pl
from src.common.dataset_catalog.api.get_datasets import get_mouse_tec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.manuscript.figures.layout.panel_layout import PanelLayout
from src.manuscript.panels.fig_02.genomic_windows.create_genomic_clusters import _plot_cluster
from src.manuscript.shared.config.image_config import DPI
from src.manuscript.shared.fig.fig import create_figure, save_svg

mouse_dataset_id: str = "mouse_gse8564_balbc_mtec_high"

show_clusters: dict[str, str] = {
    "Mast cell proteases": "Mcpt1|Mcpt2|Mcpt4|Mcpt8|Ctsg|Gzme",
    "Kallikreins I": "Klk8|Klk7|Klk6|Klk5|Klk4|Klk1b7-ps|Klk1b8|Klk1b1",
    "Kallikreins II": "Klk1b21|Klk1b16|Klk1b24|Klk1b3|Klk1b4|Klk1b5|Klk1",
    "Insulin locus": "Ins2|Th|Ascl2",
    "OAS (ISGs)": "Oas1e|Oas1c|Oas1b|Oas1f|Oas1h",
    "GBP (ISGs)": "Gbp7|Gbp3|Gbp2b|Gbp2",
}


def create_mice_window(layout: PanelLayout, output_file: Path) -> Path:
    mouse_datasets: list[DatasetMetadata] = get_mouse_tec_high()

    mouse_dataset: DatasetMetadata = next(
        dataset
        for dataset in mouse_datasets
        if dataset.dataset_id == mouse_dataset_id
    )

    selected_genes = [gene for genes in show_clusters.values() for gene in genes.split("|")]
    tec_df: pl.DataFrame = (
        mouse_dataset.get_standardized_df()
        .filter(pl.col("gene_symbol").is_in(selected_genes))
        .with_columns(pl.lit(mouse_dataset_id).alias("dataset_id"))
    )

    fig = create_figure(width_mm=layout.width, height_mm=layout.height, dpi=DPI)
    fig.set_layout_engine("constrained")
    axes = fig.subplots(2, 3).flat

    for ax, (cluster_name, genes) in zip(axes, show_clusters.items()):
        _plot_cluster(ax, cluster_name, genes.split("|"), tec_df, dataset_id=mouse_dataset_id)

    save_svg(fig, output_file)
    return output_file
