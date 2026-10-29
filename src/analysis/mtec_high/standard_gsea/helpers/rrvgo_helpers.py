from pathlib import Path
from typing import Literal, cast

from config.paths import TEMP_FOLDER
from src.analysis.mtec_high.standard_gsea.helpers.apear_helpers import get_go_abbreviation
from src.analysis.mtec_high.standard_gsea.helpers.simplify_enrichment_helpers import get_reference_species
from src.common.feratures.rrvgo.run_rrvgo import run_rrvgo, RrvgoInput
from src.common.paths import PROJECT_ROOT
import polars as pl

def run_rrvgo_direction(
    *,
    direction: str,
    gsea_mode: str,
    gsea_files: list[Path],
    output_folder: Path,
    padj_cutoff: float,
) -> None:
    temp_folder: Path = PROJECT_ROOT / TEMP_FOLDER / "rrvgo"
    temp_folder.mkdir(parents=True, exist_ok=True)

    for gsea_file in gsea_files:
        go_abbreviation: str = get_go_abbreviation(
            gsea_file,
            gsea_mode=gsea_mode,
        )
        gsea_df: pl.DataFrame = pl.read_csv(gsea_file)

        if "dataset_id" not in gsea_df.columns:
            raise ValueError(f"GSEA input is missing the dataset_id column: {gsea_file}")
        if gsea_df.get_column("dataset_id").null_count() > 0:
            raise ValueError(f"GSEA input contains a null dataset_id: {gsea_file}")

        for (dataset_id_value,), dataset_df in gsea_df.group_by(
            "dataset_id",
            maintain_order=True,
        ):
            dataset_id: str = str(dataset_id_value)
            output_prefix: str = f"{dataset_id}_{go_abbreviation}"
            dataset_gsea_file: Path = (
                temp_folder / f"{output_prefix}_{direction}.csv"
            )
            dataset_output_folder: Path = output_folder / output_prefix

            dataset_df.write_csv(dataset_gsea_file)
            dataset_output_folder.mkdir(parents=True, exist_ok=True)

            significant_go_terms = (
                dataset_df
                .filter(
                    pl.col("fdr_q_value").is_not_null()
                    & (pl.col("fdr_q_value") <= padj_cutoff)
                    & pl.col("term").str.contains(r"^\s*GO:[0-9]{7}")
                )
                .select(
                    pl.col("term")
                    .str.extract(r"^\s*(GO:[0-9]{7})", 1)
                    .alias("go_id")
                )
                .unique()
            )

            if significant_go_terms.height < 2:
                print(f"Skipping {dataset_id} {go_abbreviation} {direction}: only {significant_go_terms.height} significant GO term(s).")
                continue

            run_rrvgo(
                RrvgoInput(
                    gsea_file=dataset_gsea_file,
                    reference_species=get_reference_species(dataset_id),
                    ontology=cast(Literal["BP", "MF", "CC"], go_abbreviation),
                    padj_cutoff=padj_cutoff,
                    output_folder=dataset_output_folder,
                )
            )
