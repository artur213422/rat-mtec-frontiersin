from pathlib import Path
import polars as pl
from matplotlib import pyplot as plt

from src.common.dataset_catalog.api.get_datasets import get_mouse_tec_low
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.feratures.visualization.venn_diagram.venn_diagram import create_venn_diagram
from src.common.feratures.visualization.venn_diagram.venn_diagram_types import VennDiagramSet, VennDiagramData, \
    VennDiagramStyle, VennDiagramSpec
from src.common.paths import find_file_in_project, resolve_dir_from_root
from src.exploration.mtec_high.general_statistics.utils.prepare_venn_diagram import write_overlapped

##########################################################################################
#                                       PREPARE DATA
##########################################################################################
ortholog_map_file: Path = find_file_in_project("HOM_AllOrganism.rpt.csv")

output_folder: Path = resolve_dir_from_root("data/mtec_low/analysis/affy_vs_rna_seq")

mouse_datasets: list[DatasetMetadata] = get_mouse_tec_low()

diagrams: dict[str, list[str]] = {
    "rna_seq":
        ["mouse_gse222285_mtec_low", "mouse_gse65617_mtec_low"],
    "they should be same":
        ["mouse_gse222285_mtec_low", "mouse_gse205410_mtec_low"],
    "affy_vs_rna_seq":
        ["mouse_gse2585_mtec_low", "mouse_gse222285_mtec_low", "mouse_gse65617_mtec_low"],
}

DATASET_COLORS: dict[str, str] = {
    "mouse_gse2585_mtec_low": "#4C78A8",
    "mouse_gse222285_mtec_low": "#F58518",
    "mouse_gse205410_mtec_low": "#54A24B",
    "mouse_gse65617_mtec_low": "#E45756",
}

dataset_by_id: dict[str, DatasetMetadata] = {
    dataset.dataset_id: dataset
    for dataset in mouse_datasets
}


##########################################################################################
#                                       EXTRACT GENE SYMBOLS
##########################################################################################
gene_symbols_dataset_id: dict[str, list[str]] = {}

for dataset in mouse_datasets:
    tec_df: pl.DataFrame = dataset.get_standardized_df()

    gene_symbols_dataset_id[dataset.dataset_id] = tec_df.get_column("gene_symbol").to_list()


##########################################################################################
#                                       PROCESS FILE
##########################################################################################


for diagram_label, dataset_ids in diagrams.items():
    print(f"Creating diagram: {diagram_label}")

    venn_sets: list[VennDiagramSet] = []

    for dataset_id in dataset_ids:
        try:
            dataset_genes: list[str] = gene_symbols_dataset_id[dataset_id]
            dataset_metadata: DatasetMetadata = dataset_by_id[dataset_id]
            dataset_color: str = DATASET_COLORS[dataset_id]
        except KeyError as error:
            raise ValueError(
                f"Dataset '{dataset_id}' required by diagram "
                f"'{diagram_label}' was not found."
            ) from error

        venn_sets.append(
            VennDiagramSet(
                label=dataset_metadata.dataset_name,
                data=dataset_genes,
                color=dataset_color,
            )
        )

    data = VennDiagramData(
        sets=tuple(venn_sets),
    )

    spec = VennDiagramSpec(
        data=data,
        style=VennDiagramStyle(),
    )

    write_overlapped(
        data,
        output_folder / f"{diagram_label}.txt",
    )

    fig, ax = plt.subplots(figsize=(10.0, 10.0))

    create_venn_diagram(
        ax=ax,
        spec=spec,
    )

    fig.suptitle(diagram_label.replace("_", " "))
    fig.tight_layout()

    ##########################################################################################
    #                                       SAVE
    ##########################################################################################
    output_file: Path = output_folder / f"{diagram_label}.png"

    fig.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Saved: {output_file}")