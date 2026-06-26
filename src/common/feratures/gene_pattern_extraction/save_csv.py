from pathlib import Path

import polars as pl

from src.common.feratures.gene_pattern_extraction.models import SharedGeneReferences, GeneReferencePerDataset, \
    GeneReference


def save_gene_references(
    output_dir: Path,
    mapped_chemokines: list[GeneReferencePerDataset],
    shared_chemokines: SharedGeneReferences,
    pattern_name: str,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str | int]] = []

    for dataset in mapped_chemokines:
        for ref in dataset.gen_refs:
            rows.append(
                {
                    "dataset_id": dataset.dataset_id,
                    "gene_symbol": ref.gene_symbol,
                    "ortholog_id": ref.ortholog_id,
                }
            )

    pl.DataFrame(rows).write_csv(
        output_dir / f"{pattern_name}_per_dataset.csv"
    )

    shared_rows: list[dict[str, str | int]] = [
        {
            "gene_symbol": chemokine.gene_symbol,
            "ortholog_id": chemokine.ortholog_id,
        }
        for chemokine in shared_chemokines.shared_refs
    ]

    pl.DataFrame(shared_rows).write_csv(
        output_dir / f"shared_{pattern_name}.csv"
    )

    missing_rows: list[dict[str, str | int]] = [
        {
            "gene_symbol": chemokine.gene_symbol,
            "ortholog_id": chemokine.ortholog_id,
        }
        for chemokine in shared_chemokines.missing_refs
    ]

    pl.DataFrame(missing_rows).write_csv(
        output_dir / f"missing_{pattern_name}.csv"
    )


def _format_chemokine_list(refs: list[GeneReference]) -> str:
    if not refs:
        return "None"

    formatted: list[str] = [
        f"{ref.gene_symbol} ({ref.ortholog_id})"
        for ref in sorted(
            refs,
            key=lambda item: item.gene_symbol,
        )
    ]

    return ", ".join(formatted)


def write_report(
    report_file: Path,
    mapped_chemokines: list[GeneReferencePerDataset],
    shared: SharedGeneReferences,
    pattern_name: str,
) -> None:
    report_file.parent.mkdir(parents=True, exist_ok=True)

    with report_file.open("w", encoding="utf-8") as file:
        file.write(f"{pattern_name.upper()} PER DATASET\n")
        file.write("======================\n\n")

        for dataset in mapped_chemokines:
            file.write(
                f"{dataset.dataset_id} "
                f"(n={len(dataset.gen_refs)}):\n"
            )

            file.write(
                "  " + _format_chemokine_list(dataset.gen_refs)
            )

            file.write("\n\n")

        file.write(
            f"SHARED {pattern_name.upper()} "
            f"(n={len(shared.shared_refs)}):\n"
        )

        file.write(
            "  " + _format_chemokine_list(shared.shared_refs)
        )

        file.write("\n\n")

        file.write(
            f"MISSING {pattern_name.upper()} "
            f"(n={len(shared.missing_refs)}):\n"
        )

        file.write(
            "  " + _format_chemokine_list(shared.missing_refs)
        )

        file.write("\n")