from pathlib import Path
import polars as pl

from src.common.models.dataset_metadata.model import DatasetCatalog, DatasetMetadata
from src.common.models.dataset_metadata.ops import build_dataset_catalog
from src.common.paths import find_file_in_project, resolve_dir_from_root
from src.data_preprocessing.internal.models.ortholog_mapping_result.model import OrthologMappingResult
from src.data_preprocessing.internal.models.ortholog_mapping_result.ops import save_ortholog_mapping_result
from src.data_preprocessing.internal.orthology_mapping import find_ortholog, OrthologRequest, OrthologSpecies, Status, \
    drop_genes_with_multi_orthologs_and_save

##########################################################################################
#                                       PREPARE DATA
##########################################################################################
ortholog_map_file: Path = find_file_in_project("HOM_AllOrganism.rpt.csv")

datasets: DatasetCatalog = build_dataset_catalog()
mouse_datasets: list[DatasetMetadata] = datasets.get_mouse_datasets()

output: Path = resolve_dir_from_root("data/standardized/mouse")
source_path: Path = resolve_dir_from_root("data/intermediate/mouse")
ortholog_map: pl.DataFrame = pl.read_csv(ortholog_map_file)


##########################################################################################
#                                       PROCESS FILE
##########################################################################################

for dataset in mouse_datasets:
    ### Prepare data
    print(f"Processing {dataset.file.stem}")

    source_file: Path = source_path / f"{dataset.file.stem}.csv"

    if not source_file.exists():
        raise FileNotFoundError(f"{source_file} does not exist")

    tec_df: pl.DataFrame = pl.read_csv(source_file)

    mapping_result: OrthologMappingResult = OrthologMappingResult(dataset.dataset_id)
    mapping_result.total_genes = tec_df.height

    mapped_rows: list[dict[str, object]] = []

    ### Map orthologs
    for row in tec_df.iter_rows(named=True):
        gene: str = row["gene_symbol"]
        result: OrthologRequest = find_ortholog(gene, OrthologSpecies.MOUSE, ortholog_map)

        match result.status:
            case Status.OK:
                if len(result.rat_symbols) > 1:
                    mapping_result.rat_polymorphism_count += 1

                if len(result.mouse_symbols) > 1:
                    mapping_result.mouse_polymorphism_count += 1


            case Status.ORTHOLOGY_ID_NOT_FOUND:
                mapping_result.genes_with_no_orthology += 1

        mapped_rows.append({
            **row,

            # new fields
            "ortholog_id": result.ortholog_id,
            "rat_ortholog_name": ";".join(result.rat_symbols),

            # optional but useful audit fields
            "mouse_ortholog_names": ";".join(result.rat_symbols),
            "n_mouse_orthologs": len(result.mouse_symbols),
            "n_rat_orthologs": len(result.rat_symbols),
        })

    tec_df = pl.DataFrame(mapped_rows)

    ## Drop unmapped rows
    tec_df = tec_df.filter(
        pl.col("ortholog_id").is_not_null()
    )


    mapping_result.genes_with_multiple_orthologs = tec_df.filter(
        pl.col("ortholog_id").str.contains("///", literal=True)
    ).height

    ## Drop genes with multiple ortholog IDs
    dropped_orthologs: Path = output / f"{dataset.dataset_id}_bad_orthologs.csv"
    tec_df: pl.DataFrame = drop_genes_with_multi_orthologs_and_save(tec_df, dropped_orthologs)

    mapping_result.mapped_genes_count = tec_df.height

    ##########################################################################################
    #                                       SAVE
    ##########################################################################################

    cleaned_output_file: Path = output / f"{dataset.dataset_id}.csv"
    tec_df.write_csv(cleaned_output_file)

    save_ortholog_mapping_result(
        result=mapping_result,
        output_path=output,
        filename=dataset.dataset_id,
    )
print("Done!")