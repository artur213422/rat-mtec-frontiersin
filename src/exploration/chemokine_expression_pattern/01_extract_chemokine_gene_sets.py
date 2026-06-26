import re
from pathlib import Path
import polars as pl

from src.common.feratures.gene_pattern_extraction.methods import find_shared_gene_references, \
    find_gene_references_by_pattern
from src.common.models.dataset_metadata.model import DatasetCatalog
from src.common.models.dataset_metadata.ops import build_dataset_catalog
from src.common.feratures.gene_pattern_extraction.models import GeneReference, GeneReferencePerDataset, \
    SharedGeneReferences
from src.common.paths import resolve_dir_from_root, find_file_in_project


from src.common.feratures.gene_pattern_extraction.save_csv import write_report, save_gene_references

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


source: Path = resolve_dir_from_root("data/standardized")
output_dir: Path = resolve_dir_from_root("data/analysis/chemokine_expression_pattern/chemokine_lists")
report_file: Path = output_dir / "chemokine_coverage_report.txt"

datasets: DatasetCatalog = build_dataset_catalog()

chemokine_patterns: list[re.Pattern[str]] = [
    re.compile(r"^Ccl\d+[a-z]?$", re.IGNORECASE),
    re.compile(r"^Cxcl\d+[a-z]?$", re.IGNORECASE),
    re.compile(r"^Xcl\d+[a-z]?$", re.IGNORECASE),
    re.compile(r"^Cx3cl\d+[a-z]?$", re.IGNORECASE),
]

mapped_chemokines: list[GeneReferencePerDataset] = []
##########################################################################################
#                                       PREPARE DATA
##########################################################################################

for dataset in datasets.datasets:
    dataset_file_name = f"{dataset.dataset_id}.csv"

    file: Path = find_file_in_project(dataset_file_name, source)

    tec_df: pl.DataFrame = pl.read_csv(file)

    chemokines: list[GeneReference] = find_gene_references_by_pattern(
        data=tec_df,
        gene_patterns=chemokine_patterns,
        gene_col="gene_symbol",
        ortholog_col="ortholog_id",
    )

    mapped_chemokines.append(
        GeneReferencePerDataset(
            dataset_id=dataset.dataset_id,
            gen_refs=chemokines,
        )
    )

    print(dataset.dataset_id)
    print([hit.gene_symbol for hit in chemokines])

shared_chemokines: SharedGeneReferences = find_shared_gene_references(mapped_chemokines)
##########################################################################################
#                                       SAVE
##########################################################################################
save_gene_references(output_dir, mapped_chemokines, shared_chemokines, "chemokines")

write_report(report_file, mapped_chemokines, shared_chemokines, "chemokines")
