from pathlib import Path
from tempfile import TemporaryDirectory

import polars as pl

from src.analysis.mtec_high.isg_pathways_gsea.shared.gsea_helpers import prepare_isg_gene_sets
from src.common.dataset_catalog.api.get_datasets import get_mtec_high
from src.common.dataset_catalog.internals.types import DatasetMetadata
from src.common.feratures.padog.run_padog import PadogInput, run_padog
from src.common.paths import PROJECT_ROOT

##########################################################################################
# PREPARE DATA
##########################################################################################

input_output: dict[Path, Path] = {
    Path("resources/interferon/gse161664_lung/standardized/picked_isgs/gse161664_ifn_specific_induced_isgs.csv"):
        Path("data/mtec_high/derived_data/isg_pathways/gse161664_human_lung_epithelia_cells/padog/ifn_specific_induced_isgs"),
    Path("resources/interferon/gse161664_lung/standardized/picked_isgs/gse161664_ifn_specific_intersection_isgs.csv"):
        Path("data/mtec_high/derived_data/isg_pathways/gse161664_human_lung_epithelia_cells/padog/ifn_vs_ctrl_intersection"),
    Path("resources/interferon/gse161664_lung/standardized/picked_isgs/gse161664_ifn_vs_control_isgs.csv"):
        Path("data/mtec_high/derived_data/isg_pathways/gse161664_human_lung_epithelia_cells/padog/ifn_vs_ctrl"),
}

def write_gmt(gmt_file: Path, gene_sets: dict[str, list[str]]) -> None:
    gmt_lines: list[str] = []

    for term, genes in sorted(gene_sets.items()):
        clean_term: str = term.strip()
        clean_genes: list[str] = list(dict.fromkeys(gene.strip() for gene in genes if gene.strip()))

        if not clean_term:
            raise ValueError("Gene-set names cannot be empty")
        if any(character in clean_term for character in "\t\r\n"):
            raise ValueError(f"Invalid character in gene-set name: {term!r}")
        if any(
            any(character in gene for character in "\t\r\n")
            for gene in clean_genes
        ):
            raise ValueError(f"Invalid character in gene set: {term!r}")

        gmt_lines.append("\t".join([clean_term, "NA", *clean_genes]))

    if not gmt_lines:
        raise ValueError("No gene sets were provided")

    gmt_file.write_text("\n".join(gmt_lines) + "\n", encoding="utf-8")


def prepare_expression_matrix(dataset: DatasetMetadata) -> pl.DataFrame:
    expression_df: pl.DataFrame = (
        dataset
        .get_standardized_orthologs_df()
        .select(
            pl.col("ortholog_id").cast(pl.String).str.strip_chars(),
            *dataset.get_control_group(),
            *dataset.get_ko_group(),
        )
        .drop_nulls()
        .filter(pl.col("ortholog_id") != "")
    )

    duplicate_ids: pl.DataFrame = (
        expression_df
        .group_by("ortholog_id")
        .len()
        .filter(pl.col("len") > 1)
    )
    if duplicate_ids.height:
        raise ValueError(f"{dataset.dataset_id}: expression matrix contains duplicate ortholog IDs. Examples:\n{duplicate_ids.head(20)}")

    return expression_df



datasets: list[DatasetMetadata] = get_mtec_high()
temp_root: Path = PROJECT_ROOT / "data/temp/padog"
temp_root.mkdir(parents=True, exist_ok=True)

######################################################################################
# PROCESS DATA
######################################################################################

for relative_input_file, relative_out_dir in input_output.items():
    input_file: Path = PROJECT_ROOT / relative_input_file
    output_file: Path = PROJECT_ROOT / relative_out_dir / "padog_results.csv"
    output_file.parent.mkdir(parents=True, exist_ok=True)

    gene_sets: dict[str, list[str]] = prepare_isg_gene_sets(input_file)

    if len(gene_sets) < 3:
        reason = (
            f"PADOG unavailable: requires at least 3 gene sets, "
            f"but {input_file.name} contains {len(gene_sets)}."
        )
        output_file.with_name("padog_status.txt").write_text(
            reason + "\n", encoding="utf-8"
        )
        print(reason)
        continue

    experiment_results: list[pl.DataFrame] = []

    # All PADOG inputs and per-dataset outputs are removed on success or error.
    with TemporaryDirectory(dir=temp_root,prefix=f"{relative_out_dir.name}_",) as temporary_directory:
        temp_dir = Path(temporary_directory)
        gmt_file: Path = temp_dir / "isg_gene_sets.gmt"
        write_gmt(gmt_file, gene_sets)

        for dataset in datasets:
            sample_count: int = (len(dataset.get_control_group()) + len(dataset.get_ko_group()))
            if sample_count <= 4:
                print(f"Skipping {dataset.dataset_id}: PADOG requires at least 5 samples")
                continue

            print("Analysing dataset:", dataset.dataset_id)

            expression_matrix_file: Path = temp_dir / f"{dataset.dataset_id}_expression.csv"

            dataset_output_file: Path = temp_dir / f"{dataset.dataset_id}_padog.csv"


            prepare_expression_matrix(dataset).write_csv(expression_matrix_file)

            run_padog(
                PadogInput(
                    expression_matrix_file=expression_matrix_file,
                    gmt_file=gmt_file,
                    output_file=dataset_output_file,
                    control_group=dataset.get_control_group(),
                    ko_group=dataset.get_ko_group(),
                )
            )

            if not dataset_output_file.is_file():
                raise RuntimeError(f"PADOG did not create a result for {dataset.dataset_id}: {dataset_output_file}")

            experiment_results.append(
                pl.read_csv(dataset_output_file).with_columns(
                    pl.lit(dataset.dataset_id).alias("dataset_id")
                )
            )

    ##################################################################################
    # SAVE
    ##################################################################################

    if not experiment_results:
        raise RuntimeError(f"PADOG produced no results for {input_file}")

    combined_results: pl.DataFrame = pl.concat(
        experiment_results,
        how="vertical",
    ).select("dataset_id", pl.exclude("dataset_id"))

    combined_results.write_csv(output_file)
    print(f"Saved combined PADOG results: {output_file}")



