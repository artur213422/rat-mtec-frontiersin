from pathlib import Path

from matplotlib.axes import Axes
from matplotlib.figure import Figure

from src.common.dataset_catalog.api.get_datasets import read_file_with_ortholog_id
from src.common.paths import find_file_in_project, PROJECT_ROOT
from src.manuscript.figures.layout.panel_layout import PanelLayout
from src.manuscript.shared.config.dataset_labels import DATASET_LABELS
from src.manuscript.shared.config.image_config import ManuscriptColors, DPI
from src.manuscript.shared.fig.fig import create_figure, save_svg
import polars as pl

from src.manuscript.shared.heatmap.draw_heatmap import draw_manuscript_heatmap
from src.manuscript.shared.heatmap.heatmap_input import ManuscriptHeatmapInput
from src.resource_processing.standard_gsea_pathways.internal.enums import Namespace

gsea_folder: Path = PROJECT_ROOT / "data/mtec_high/derived_data/standard_gsea/both_sided_gsea"
gsea_files = {
    Namespace.BP: find_file_in_project("go_biological_process_gsea.csv", gsea_folder),
    Namespace.MF: find_file_in_project("go_molecular_function_gsea.csv", gsea_folder),
    Namespace.CP: find_file_in_project("go_cellular_component_gsea.csv", gsea_folder),

}

PATHWAYS = [
    (Namespace.BP, "GO:0034340", "Response to type I interferon"),
    (Namespace.BP, "GO:0034341", "Response to type II interferon"),

    (Namespace.BP, "GO:0051607", "Defense response to virus"),
    (Namespace.BP, "GO:0019730", "Antimicrobial humoral response"),
    (Namespace.BP, "GO:0002495", "MHC-II peptide antigen processing and presentation"),

    (Namespace.BP, "GO:0006260", "DNA replication"),
    (Namespace.BP, "GO:0042254", "Ribosome biogenesis"),
    (Namespace.BP, "GO:0051147", "Regulation of muscle cell differentiation"),


    (Namespace.BP, "GO:0031424", "keratinization"),
    (Namespace.BP, "GO:0043588", "Skin development"),
    (Namespace.BP, "GO:0042574", "Retinal metabolic process"),

    (Namespace.BP, "GO:0007218", "Neuropeptide signaling pathway"),
    (Namespace.BP, "GO:0046883", "Regulation of hormone secretion"),
    (Namespace.BP, "GO:0051384", "Response to glucocorticoid"),
    (Namespace.MF, "GO:0004252", "Serine-type endopeptidase activity"),
]

NES_LIMIT = 3.0

def _prepare_gsea_data() -> pl.DataFrame:
    pathway_labels: dict[str, str] = {
        go_id: label for _, go_id, label in PATHWAYS
    }

    if not PATHWAYS or len(pathway_labels) != len(PATHWAYS):
        raise ValueError("PATHWAYS must contain unique GO IDs.")

    if len(set(pathway_labels.values())) != len(PATHWAYS):
        raise ValueError("PATHWAYS must contain unique display labels.")

    required: set[str] = {
        "term",
        "dataset_id",
        "nes",
        "nominal_p_value",
        "fdr_q_value",
    }

    namespace_data: list[pl.DataFrame] = []

    for namespace in dict.fromkeys(namespace for namespace, _, _ in PATHWAYS):
        if namespace not in gsea_files:
            raise ValueError(f"Unsupported GO namespace: {namespace!r}")

        go_gsea: pl.DataFrame = read_file_with_ortholog_id(gsea_files[namespace])
        missing_columns: set[str] = required - set(go_gsea.columns)

        if missing_columns:
            raise ValueError(
                f"Missing GSEA columns for {namespace.value}: {sorted(missing_columns)}"
            )

        namespace_go_ids: list[str] = [
            go_id for term_namespace, go_id, _ in PATHWAYS
            if term_namespace == namespace
        ]

        namespace_data.append(
            go_gsea
            .with_columns(
                pl.col("term")
                .str.split("|")
                .list.first()
                .str.strip_chars()
                .alias("go_id")
            )
            .filter(
                pl.col("go_id").is_in(namespace_go_ids),
                pl.col("dataset_id").is_in(list(DATASET_LABELS)),
            )
            .select(
                pl.col("go_id"),
                pl.col("dataset_id").cast(pl.String),
                pl.col("nes").cast(pl.Float64, strict=False),
                pl.col("nominal_p_value").cast(pl.Float64, strict=False),
                pl.col("fdr_q_value").cast(pl.Float64, strict=False),
            )
        )

    data: pl.DataFrame = pl.concat(namespace_data)

    missing_pathways: set[str] = set(pathway_labels) - set(
        data.get_column("go_id").to_list()
    )

    if missing_pathways:
        raise ValueError(
            f"GO pathways not found in GSEA data: {sorted(missing_pathways)}"
        )

    duplicates: pl.DataFrame = (
        data
        .group_by("go_id", "dataset_id")
        .len()
        .filter(pl.col("len") > 1)
    )

    if not duplicates.is_empty():
        raise ValueError(f"Duplicate GO/dataset results:\n{duplicates}")

    return (
        data
        .with_columns(
            pl.col("go_id")
            .replace_strict(pathway_labels)
            .alias("pathway")
        )
        .select(
            "pathway",
            "dataset_id",
            "nes",
            "nominal_p_value",
            "fdr_q_value",
        )
    )


def _draw_gsea_heatmap(ax: Axes, data: pl.DataFrame) -> None:
    dataset_ids: list[str] = [
        dataset_id
        for dataset_id in DATASET_LABELS
        if dataset_id in data.get_column("dataset_id").to_list()
    ]

    manuscript_input: ManuscriptHeatmapInput = ManuscriptHeatmapInput(
        data=data,
        value_col="nes",
        row_label_col="pathway",
        dataset_id_col="dataset_id",
        p_col="nominal_p_value",
        fdr_col="fdr_q_value",

        col_order=dataset_ids,
        col_labels=DATASET_LABELS,
        row_order=[label for _, _, label in PATHWAYS],

        italic_row_labels=False,
        show_values=True,

        vmin=-NES_LIMIT,
        vmax=NES_LIMIT,
        min_color=ManuscriptColors.BLUE,
        max_color=ManuscriptColors.RED,
        na_color=ManuscriptColors.WHITE,

        p_threshold=0.05,
        fdr_threshold=0.05,

        colorbar_label="NES"
    )

    draw_manuscript_heatmap(ax, manuscript_input)

def create_gsea_table(layout: PanelLayout, output_file: Path) -> Path:
    fig: Figure = create_figure(width_mm=layout.width, height_mm=layout.height, dpi=DPI)

    ax = fig.add_axes((0, 0, 1, 1))

    data: pl.DataFrame = _prepare_gsea_data()
    _draw_gsea_heatmap(ax, data)

    save_svg(fig, output_file)

    return output_file
