from pathlib import Path

import polars as pl

from src.analysis.mtec_high.genomics.curated_loci._02_map_curated_loci.gene_family.helpers.select_gene_family_members import select_gene_family_members
from src.analysis.mtec_high.genomics.curated_loci.curated_clusters import GENE_FAMILIES, CuratedLocus
from src.analysis.mtec_high.genomics.curated_loci.utils.configs.paths import CURATED_LOCI_FILE
from src.analysis.mtec_high.genomics.curated_loci.utils.configs.types import LocusDefinitionType
from src.analysis.mtec_high.genomics.curated_loci.utils.parsers.handle_curated_loci_file import build_curated_locus, overwrite_to_curated_loci_file
from src.common.dataset_catalog.api.get_datasets import read_file_with_ortholog_id
from src.common.paths import find_file_in_project

# ----------------------------------------------------------------------------------------
# Prepare Data
# ----------------------------------------------------------------------------------------

clusters_file: Path = find_file_in_project("datasets_with_clusters.csv")
clusters_df: pl.DataFrame = read_file_with_ortholog_id(clusters_file)

dataset_frames: list[pl.DataFrame] = clusters_df.partition_by("dataset_id", maintain_order=True)

# ----------------------------------------------------------------------------------------
# Pick Families
# ----------------------------------------------------------------------------------------

loci: list[pl.DataFrame] = []

for locus_name, gene_family in GENE_FAMILIES.items():
    gene_family: CuratedLocus

    print(f"Mapping curated locus: {locus_name}")

    for tec_df in dataset_frames:
        dataset_id: str = tec_df.get_column("dataset_id").item(0)

        member_df: pl.DataFrame = select_gene_family_members(tec_df, gene_family.gene_patterns)

        if member_df.is_empty():
            print(f"  No {locus_name} members found in {dataset_id}")
            continue

        print(f"  {dataset_id}: {member_df.height} family members found")

        locus_df: pl.DataFrame = build_curated_locus(
            tec_df=tec_df,
            member_df=member_df,
            locus_name=locus_name,
            locus_type=gene_family.definition_type,
            dataset_id=dataset_id,
        )

        loci.append(locus_df)

# ----------------------------------------------------------------------------------------
# Save
# ----------------------------------------------------------------------------------------

if not loci:
    raise ValueError("No curated gene-family members were found. Check GENE_FAMILIES and the input gene_symbol values.")

result: pl.DataFrame = pl.concat(loci, how="vertical")

overwrite_to_curated_loci_file(
    result=result,
    locus_definition_type=LocusDefinitionType.GENE_FAMILY,
)

print(f"Saved {result.height} gene-family rows across {result['dataset_id'].n_unique()} datasets to:\n {CURATED_LOCI_FILE.resolve()}")