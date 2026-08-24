from pathlib import Path

import polars as pl
from matplotlib import pyplot as plt

from src.common.dataset_catalog.api.get_datasets import (
    get_mouse_tec_low,
    get_rat_tec_low,
)
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.feratures.visualization.venn_diagram.venn_diagram import (
    create_venn_diagram,
)
from src.common.feratures.visualization.venn_diagram.venn_diagram_types import (
    VennDiagramData,
    VennDiagramSet,
    VennDiagramSpec,
    VennDiagramStyle,
)
from src.common.paths import resolve_dir_from_root
from src.common.types.enums import Colors
from src.exploration.mtec_high.general_statistics.utils.prepare_venn_diagram import (
    write_overlapped,
)


##########################################################################################
#                                       PREPARE DATA
##########################################################################################

output_folder: Path = resolve_dir_from_root(
    "data/mtec_low/analysis/affy_vs_rna_seq/venn_diagram/rat_overlap"
)

output_folder.mkdir(
    parents=True,
    exist_ok=True,
)

mouse_datasets: list[DatasetMetadata] = get_mouse_tec_low()
rat_dataset: DatasetMetadata = get_rat_tec_low()

all_datasets: list[DatasetMetadata] = [
    *mouse_datasets,
    rat_dataset,
]

diagrams: dict[str, list[str]] = {
    "rna_seq": [
        "mouse_gse222285_mtec_low",
        "mouse_gse65617_mtec_low",
    ],
    "rat_vs_affy_mouse": [
        "rat_mtec_low",
        "mouse_gse2585_mtec_low",
    ],
    "rat_vs_rna_seq": [
        "rat_mtec_low",
        "mouse_gse222285_mtec_low",
        "mouse_gse65617_mtec_low",
    ],
}

DATASET_COLORS: dict[str, str] = {
    "mouse_gse2585_mtec_low": "#4C78A8",
    "mouse_gse222285_mtec_low": "#F58518",
    "mouse_gse65617_mtec_low": "#E45756",
    "rat_mtec_low": (
        Colors.RED.value
        if isinstance(Colors.RED, Colors)
        else str(Colors.RED)
    ),
}


##########################################################################################
#                              INDEX DATASETS AND GENES
##########################################################################################

dataset_by_id: dict[str, DatasetMetadata] = {
    dataset.dataset_id: dataset
    for dataset in all_datasets
}

gene_symbols_by_dataset_id: dict[str, list[str]] = {}

for dataset in all_datasets:
    tec_df: pl.DataFrame = dataset.get_standardized_df()

    if "gene_symbol" not in tec_df.columns:
        raise ValueError(
            f"Dataset '{dataset.dataset_id}' does not contain "
            "a 'gene_symbol' column."
        )

    gene_symbols: list[str] = (
        tec_df
        .select(
            pl.col("gene_symbol")
            .cast(pl.String)
            .str.strip_chars()
        )
        .filter(
            pl.col("gene_symbol").is_not_null()
            & (pl.col("gene_symbol") != "")
        )
        .unique()
        .sort("gene_symbol")
        .get_column("gene_symbol")
        .to_list()
    )

    gene_symbols_by_dataset_id[
        dataset.dataset_id
    ] = gene_symbols

    print(
        f"{dataset.dataset_id}: "
        f"{len(gene_symbols):,} unique gene symbols"
    )


##########################################################################################
#                                   VALIDATE SETTINGS
##########################################################################################

required_dataset_ids: set[str] = {
    dataset_id
    for dataset_ids in diagrams.values()
    for dataset_id in dataset_ids
}

missing_metadata: set[str] = (
    required_dataset_ids
    - set(dataset_by_id)
)

missing_gene_data: set[str] = (
    required_dataset_ids
    - set(gene_symbols_by_dataset_id)
)

missing_colors: set[str] = (
    required_dataset_ids
    - set(DATASET_COLORS)
)

if missing_metadata:
    raise ValueError(
        "Missing dataset metadata for: "
        f"{sorted(missing_metadata)}"
    )

if missing_gene_data:
    raise ValueError(
        "Missing gene-symbol data for: "
        f"{sorted(missing_gene_data)}"
    )

if missing_colors:
    raise ValueError(
        "Missing dataset colors for: "
        f"{sorted(missing_colors)}"
    )


##########################################################################################
#                                  CREATE DIAGRAMS
##########################################################################################

for diagram_label, dataset_ids in diagrams.items():
    print(f"Creating diagram: {diagram_label}")

    venn_sets: list[VennDiagramSet] = []

    for dataset_id in dataset_ids:
        dataset_metadata: DatasetMetadata = (
            dataset_by_id[dataset_id]
        )

        dataset_genes: list[str] = (
            gene_symbols_by_dataset_id[dataset_id]
        )

        dataset_color: str = (
            DATASET_COLORS[dataset_id]
        )

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

    text_output_file: Path = (
        output_folder
        / f"{diagram_label}.txt"
    )

    write_overlapped(
        data,
        text_output_file,
    )

    fig, ax = plt.subplots(
        figsize=(10.0, 10.0)
    )

    create_venn_diagram(
        ax=ax,
        spec=spec,
    )

    fig.suptitle(
        diagram_label.replace("_", " ")
    )

    fig.tight_layout()

    image_output_file: Path = (
        output_folder
        / f"{diagram_label}.png"
    )

    fig.savefig(
        image_output_file,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Saved: {image_output_file}")