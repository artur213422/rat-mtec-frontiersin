import re
from pathlib import Path
import polars as pl

from src.common.feratures.gene_pattern_extraction.methods import find_gene_references_by_pattern, \
    find_shared_gene_references
from src.common.models.dataset_metadata.model import DatasetCatalog
from src.common.models.dataset_metadata.ops import build_dataset_catalog
from src.common.feratures.gene_pattern_extraction.models import GeneReference, SharedGeneReferences, \
    GeneReferencePerDataset
from src.common.paths import resolve_dir_from_root, find_file_in_project


from src.common.feratures.gene_pattern_extraction.save_csv import write_report, save_gene_references

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


source: Path = resolve_dir_from_root("data/standardized")
output_dir: Path = resolve_dir_from_root("data/analysis/kallikrein_expression_pattern/kallikrein_lists")
report_file: Path = output_dir / "kallikrein_coverage_report.txt"

datasets: DatasetCatalog = build_dataset_catalog()

kallikrein_patterns: list[re.Pattern[str]] = [
    # Canonical tissue kallikrein / kallikrein-related peptidases:
    # Klk1, Klk2, ..., Klk15
    re.compile(r"^Klk(?:[1-9]|1[0-5])$", re.IGNORECASE),

    # Rodent kallikrein 1-related peptidase family:
    # Klk1b3, Klk1b4, Klk1b22, Klk1b27, etc.
    re.compile(r"^Klk1b\d+[a-z]?$", re.IGNORECASE),

    # Plasma kallikrein:
    # Klkb1
    re.compile(r"^Klkb1$", re.IGNORECASE),
]


mapped_kallikreins: list[GeneReferencePerDataset] = []

##########################################################################################
#                                       PREPARE DATA
##########################################################################################

for dataset in datasets.datasets:
    dataset_file_name = f"{dataset.dataset_id}.csv"

    file: Path = find_file_in_project(dataset_file_name, source)

    tec_df: pl.DataFrame = pl.read_csv(file)

    chemokines: list[GeneReference] = find_gene_references_by_pattern(
        data=tec_df,
        gene_patterns=kallikrein_patterns,
        gene_col="gene_symbol",
        ortholog_col="ortholog_id",
    )

    mapped_kallikreins.append(
        GeneReferencePerDataset(
            dataset_id=dataset.dataset_id,
            gen_refs=chemokines,
        )
    )

    print(dataset.dataset_id)
    print([hit.gene_symbol for hit in chemokines])

shared_chemokines: SharedGeneReferences = find_shared_gene_references(mapped_kallikreins)

##########################################################################################
#                                       SAVE
##########################################################################################
save_gene_references(output_dir, mapped_kallikreins, shared_chemokines, "kallikreins")

write_report(report_file, mapped_kallikreins, shared_chemokines, "kallikreins")
