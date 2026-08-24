from pathlib import Path

import polars as pl
from matplotlib import pyplot as plt

from src.common.dataset_catalog.api.get_datasets import get_mtec_high, get_mtec_low
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.feratures.analysers.deg_extractor import DatasetDegs, get_degs
from src.common.feratures.visualization.venn_diagram.venn_diagram import  create_venn_diagram
from src.common.feratures.visualization.venn_diagram.venn_diagram_types import VennDiagramData, VennDiagramSet, VennDiagramSpec, VennDiagramStyle
from src.common.paths import resolve_dir_from_root
from src.common.types.enums import Colors
from src.exploration.mtec_high.general_statistics.utils.colors import get_dataset_color
from src.exploration.mtec_high.general_statistics.utils.prepare_venn_diagram import write_overlapped

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

output_folder: Path = resolve_dir_from_root("data/mtec_low/analysis/general_statistics")

datasets: list[DatasetMetadata] = get_mtec_low()

diagrams: dict[str, list[str]] = {
    "rat vs affy mouse":
        ["rat_mtec_low", "mouse_gse2585_mtec_low"],
    "rat vs rna seq":
        ["rat_mtec_low", "mouse_gse65617_mtec_low", "mouse_gse205410_mtec_low",],
    "mices":
        ["mouse_gse205410_mtec_low", "mouse_gse65617_mtec_low", "mouse_gse2585_mtec_low"],
}

##########################################################################################
#                                       EXTRACT DEGS
##########################################################################################

all_degs: list[DatasetDegs] = []

for dataset in datasets:
    print(f"Running dataset: {dataset.dataset_id}")

    tec_df: pl.DataFrame = dataset.get_standardized_df()
    color: Colors = get_dataset_color(dataset.dataset_id)

    dataset_degs: DatasetDegs = get_degs(dataset_id=dataset.dataset_id, tec_df=tec_df, color=color,)

    all_degs.append(dataset_degs)

degs_by_dataset_id: dict[str, DatasetDegs] = {
    dataset_degs.dataset_id: dataset_degs
    for dataset_degs in all_degs
}


##########################################################################################
#                                   CREATE DIAGRAMS
##########################################################################################


for diagram_label, dataset_ids in diagrams.items():
    print(f"Creating diagram: {diagram_label}")

    venn_sets: list[VennDiagramSet] = []

    for dataset_id in dataset_ids:
        try:
            dataset_degs = degs_by_dataset_id[dataset_id]
        except KeyError as error:
            raise ValueError(f"Dataset '{dataset_id}' required by diagram '{diagram_label}' was not found.") from error

        all_dataset_degs: frozenset[str] = frozenset(dataset_degs.upregulated) | frozenset(dataset_degs.downregulated)

        venn_sets.append(
            VennDiagramSet(
                label=dataset_degs.dataset_id,
                data=all_dataset_degs,
                color=dataset_degs.color,
            )
        )

    data = VennDiagramData(
            sets=tuple(venn_sets),
        )
    spec = VennDiagramSpec(
        data=VennDiagramData(
            sets=tuple(venn_sets),
        ),
        style=VennDiagramStyle(),
    )

    write_overlapped(data, output_folder / f"{diagram_label}.txt")

    fig, ax = plt.subplots(figsize=(10.0, 10.0))

    create_venn_diagram(ax=ax, spec=spec,)

    fig.suptitle(diagram_label.replace("_", " "))

    fig.tight_layout()

    output_file: Path = output_folder / f"{diagram_label}.png"

    fig.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Saved: {output_file}")