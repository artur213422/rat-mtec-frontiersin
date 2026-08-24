from pathlib import Path
import polars as pl

from src.common.models.dataset_metadata.model import DatasetCatalog, DatasetMetadata
from src.common.models.dataset_metadata.ops import build_dataset_catalog
from src.common.paths import find_file_in_project, resolve_dir_from_root
from src.data_preprocessing.mtec_high.internal.models.ortholog_mapping_result.model import OrthologMappingResult
from src.data_preprocessing.mtec_high.internal.models.ortholog_mapping_result.ops import save_ortholog_mapping_result
from src.data_preprocessing.mtec_high.internal.orthology_mapping import find_ortholog, OrthologRequest, OrthologSpecies, Status

##########################################################################################
#                                       PREPARE DATA
##########################################################################################
ortholog_map_file: Path = find_file_in_project("HOM_AllOrganism.rpt.csv")

datasets: DatasetCatalog = build_dataset_catalog()
mouse_datasets: list[DatasetMetadata] = datasets.get_mouse_datasets()

output: Path = resolve_dir_from_root("data/standardized/mouse")
source_path: Path = resolve_dir_from_root("data/raw/intermediate/cleaned_names_expression_matrix/mouse")
ortholog_map: pl.DataFrame = pl.read_csv(ortholog_map_file)


##########################################################################################
#                                       PROCESS FILE
##########################################################################################

for dataset in mouse_datasets:
    ### Prepare data
    print(f"Processing {dataset.dataset_id}")

    source_file: Path = source_path / f"{dataset.dataset_id}.csv"

    if not source_file.exists():
        raise FileNotFoundError(f"{source_file} does not exist")

    tec_df: pl.DataFrame = pl.read_csv(source_file)

    mapping_result: OrthologMappingResult = OrthologMappingResult(dataset.dataset_id)
    mapping_result.total_genes = tec_df.height

    mapped_rows: list[dict[str, object]] = []

    ### Map orthologs
    for row in tec_df.iter_rows(named=True):
        gene: str = row["gene_symbol"]
        result: OrthologRequest = find_ortholog(gene, OrthologSpecies.RAT, ortholog_map)

        has_clean_ortholog: bool = True
        match result.status:
            case Status.OK:
                if len(result.rat_symbols) > 1:
                    mapping_result.rat_polymorphism_count += 1
                    has_clean_ortholog = False
                if len(result.mouse_symbols) > 1:
                    mapping_result.mouse_polymorphism_count += 1
                    has_clean_ortholog = False

            case Status.ORTHOLOGY_ID_NOT_FOUND:
                mapping_result.genes_with_no_orthology += 1
                has_clean_ortholog = False

        mapped_rows.append({
            **row,

            # new fields
            "has_clean_ortholog": has_clean_ortholog,
            "ortholog_id": result.ortholog_id,
            "mouse_ortholog_name": ";".join(result.mouse_symbols),

            # optional but useful audit fields
            "rat_ortholog_names": ";".join(result.rat_symbols),
            "n_mouse_orthologs": len(result.mouse_symbols),
            "n_rat_orthologs": len(result.rat_symbols),
        })

    tec_df = pl.DataFrame(mapped_rows)

    mapping_result.genes_with_clean_ortholog_count = tec_df.filter(
        pl.col("has_clean_ortholog")
    ).height


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