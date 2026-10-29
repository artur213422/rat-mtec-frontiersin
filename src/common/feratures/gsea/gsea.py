import math
from pathlib import Path

from gseapy import prerank, Prerank
import polars as pl
from config.paths import RSCRIPT_EXE
from src.common.feratures.gsea.models import GseaInput, GseaResult
from src.common.paths import find_file_in_project, PROJECT_ROOT
from src.data_preprocessing.mtec_high.internal.R.run_R import run_rscript

def run_gsea(input: GseaInput) -> list[GseaResult]:
    r_script: Path = find_file_in_project("gsea.R")
    temporary_folder: Path = PROJECT_ROOT / "data/temp/gsea"

    rnk_file: Path = temporary_folder / "rnk.csv"
    gene_sets_file: Path = temporary_folder / "gene_sets.gmt"
    result_file: Path = temporary_folder / "result.csv"

    if not RSCRIPT_EXE.is_file():
        raise FileNotFoundError(f"Rscript executable was not found: {RSCRIPT_EXE}")

    if not r_script.is_file():
        raise FileNotFoundError(f"fgsea R script was not found: {r_script}")

    temporary_folder.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Write ranked genes
    # ------------------------------------------------------------------

    gene_column, score_column = input.rnk.columns

    rnk_df: pl.DataFrame = input.rnk.select(
        pl.col(gene_column)
        .str.strip_chars()
        .alias("gene"),

        pl.col(score_column)
        .cast(pl.Float64)
        .alias("score"),
    )

    invalid_rows: int = rnk_df.filter(
        pl.col("gene").is_null()
        | (pl.col("gene") == "")
        | pl.col("score").is_null()
        | pl.col("score").is_nan()
        | ~pl.col("score").is_finite()
    ).height

    if invalid_rows > 0:
        raise ValueError(f"Ranked list contains {invalid_rows} invalid rows")

    duplicate_genes: pl.DataFrame = (
        rnk_df
        .group_by("gene")
        .len()
        .filter(pl.col("len") > 1)
    )

    if duplicate_genes.height > 0:
        raise ValueError("Ranked list contains duplicate genes:\n {duplicate_genes.head(20)}")

    rnk_df.write_csv(rnk_file)

    # ------------------------------------------------------------------
    # Write gene sets in GMT format
    # ------------------------------------------------------------------

    if not input.gene_sets_dict:
        raise ValueError("No gene sets were provided")

    gmt_lines: list[str] = []

    for term, genes in sorted(input.gene_sets_dict.items()):
        clean_term: str = term.strip()

        clean_genes: list[str] = list(
            dict.fromkeys(
                gene.strip()
                for gene in genes
                if gene.strip()
            )
        )

        if not clean_term:
            raise ValueError("Gene-set names cannot be empty")

        if any(character in clean_term for character in "\t\r\n"):
            raise ValueError(f"Invalid character in gene-set name: {term!r}")

        if any(any(character in gene for character in "\t\r\n")
            for gene in clean_genes
        ):
            raise ValueError(f"Invalid character in gene set: {term!r}")

        # GMT: term, description, gene1, gene2, ...
        gmt_lines.append(
            "\t".join([clean_term, "na", *clean_genes])
        )

    gene_sets_file.write_text(
        "\n".join(gmt_lines) + "\n",
        encoding="utf-8",
    )

    # Prevent a failed R run from leaving an old result that appears valid.
    result_file.unlink(missing_ok=True)

    # ------------------------------------------------------------------
    # Run fgsea
    # ------------------------------------------------------------------

    config: dict[str, object] = {
        "min_size": input.min_size,
        "max_size": input.max_size,
        "permutation_num": input.permutation_num,
        "gene_sets_file": str(gene_sets_file),
        "result_file": str(result_file),
        "rnk_file": str(rnk_file),
        "mode": input.gsea_mode.value
    }

    run_rscript(rscript_exe=RSCRIPT_EXE, r_script=r_script, config=config)

    if not result_file.is_file():
        raise RuntimeError("fgsea did not create the expected result file: {result_file}")

    # ------------------------------------------------------------------
    # Read and validate R output
    # ------------------------------------------------------------------

    raw: pl.DataFrame = pl.read_csv(
        result_file,
        schema_overrides={
            "term": pl.String,
            "es": pl.Float64,
            "nes": pl.Float64,
            "nominal_p_value": pl.Float64,
            "fdr_q_value": pl.Float64,
            "leading_edge_genes": pl.String,
        },
    )

    required_columns: set[str] = {
        "term",
        "es",
        "nes",
        "nominal_p_value",
        "fdr_q_value",
        "leading_edge_genes",
    }

    missing_columns: set[str] = required_columns - set(raw.columns)

    if missing_columns:
        raise ValueError("fgsea result is missing columns: {sorted(missing_columns)}")

    duplicate_terms: pl.DataFrame = (
        raw
        .group_by("term")
        .len()
        .filter(pl.col("len") > 1)
    )

    if duplicate_terms.height > 0:
        raise ValueError(f"fgsea returned duplicate terms:\n{duplicate_terms}")

    results: list[GseaResult] = []

    for row in raw.iter_rows(named=True):
        term = row["term"]
        es = row["es"]
        nes = row["nes"]
        nominal_p_value = row["nominal_p_value"]
        fdr_q_value = row["fdr_q_value"]

        if not isinstance(term, str):
            raise ValueError("fgsea returned an invalid term")

        if term not in input.gene_sets_dict:
            raise ValueError(f"fgsea returned an unknown term: {term!r}")

        numeric_values = (
            es,
            nes,
            nominal_p_value,
            fdr_q_value,
        )

        if any(value is None or not math.isfinite(value) for value in numeric_values):
            raise ValueError(f"fgsea returned missing/non-finite values for {term!r}")

        if not 0 <= nominal_p_value <= 1:
            raise ValueError(f"Invalid nominal p-value for {term!r}")

        if not 0 <= fdr_q_value <= 1:
            raise ValueError(f"Invalid FDR q-value for {term!r}")

        leading_edge_value = row["leading_edge_genes"]

        leading_edge_genes: list[str] = (
            leading_edge_value.split(";")
            if isinstance(leading_edge_value, str)
            and leading_edge_value
            else []
        )

        results.append(
            GseaResult(
                term=term,
                es=float(es),
                nes=float(nes),
                nominal_p_value=float(nominal_p_value),
                fdr_q_value=float(fdr_q_value),
                leading_edge_genes=leading_edge_genes,
            )
        )

    return results

