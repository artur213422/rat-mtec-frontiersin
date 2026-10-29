from pathlib import Path

import polars as pl

from src.analysis.mtec_high.genomics.curated_loci.curated_clusters import GENE_SETS
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
# Collect Gene Set Clusters
# ----------------------------------------------------------------------------------------

loci: list[pl.DataFrame] = []

for locus_name, genes in GENE_SETS.items():
    print(f"Mapping curated locus: {locus_name}")

    for tec_df in dataset_frames:
        dataset_id: str = tec_df.get_column("dataset_id").item(0)

        member_df: pl.DataFrame = tec_df.filter(pl.col("gene_symbol").is_in(genes))

        if member_df.is_empty():
            print(f"  No {locus_name} members found in {dataset_id}")
            continue

        print(f"  {dataset_id}: {member_df.height} gene set members found")

        # Select the union of member clusters within this dataset.
        locus_df: pl.DataFrame = build_curated_locus(
            tec_df=tec_df,
            member_df=member_df,
            locus_name=locus_name,
            locus_type=LocusDefinitionType.GENE_SET,
            dataset_id=dataset_id,
        )

        loci.append(locus_df)

# ----------------------------------------------------------------------------------------
# Save
# ----------------------------------------------------------------------------------------

if loci:
    result: pl.DataFrame = pl.concat(loci, how="vertical")

    overwrite_to_curated_loci_file(result=result, locus_definition_type=LocusDefinitionType.GENE_SET)
