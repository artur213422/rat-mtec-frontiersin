from pathlib import Path

import polars as pl


def write_report(
    chemokines_df: pl.DataFrame,
    report_file: Path,
    pattern: str,
) -> None:
    required_columns: set[str] = {
        "dataset_id",
        "gene_symbol",
        "has_clean_ortholog",
    }

    missing_columns: set[str] = (
        required_columns - set(chemokines_df.columns)
    )

    if missing_columns:
        raise ValueError(
            f"Missing required columns: {sorted(missing_columns)}"
        )

    dataset_dfs: list[pl.DataFrame] = chemokines_df.partition_by(
        "dataset_id",
        maintain_order=True,
    )

    genes_per_dataset: dict[str, set[str]] = {}

    for dataset_df in dataset_dfs:
        dataset_id: str = str(
            dataset_df.get_column("dataset_id")[0]
        )

        genes: set[str] = set(
            dataset_df
            .get_column("gene_symbol")
            .drop_nulls()
            .cast(pl.String)
            .to_list()
        )

        genes_per_dataset[dataset_id] = genes

    if not genes_per_dataset:
        raise ValueError("No datasets found")

    all_gene_sets: list[set[str]] = list(
        genes_per_dataset.values()
    )

    all_genes: set[str] = set.union(*all_gene_sets)

    shared_genes: set[str] = set.intersection(
        *all_gene_sets
    )

    gene_dataset_counts: dict[str, int] = {
        gene: sum(
            gene in genes
            for genes in genes_per_dataset.values()
        )
        for gene in all_genes
    }

    unique_genes: set[str] = {
        gene
        for gene, dataset_count in gene_dataset_counts.items()
        if dataset_count == 1
    }

    report_lines: list[str] = [
        f"{pattern.capitalize()} coverage report",
        "",
        f"Datasets: {len(genes_per_dataset)}",
        f"Total unique {pattern} genes: {len(all_genes)}",
        f"Shared {pattern} genes: {len(shared_genes)}",
        f"Dataset-specific {pattern} genes: {len(unique_genes)}",
        "",
        f"Shared {pattern} genes:",
        ", ".join(sorted(shared_genes)) or "None",
        "",
        f"Dataset-specific {pattern} genes:",
        ", ".join(sorted(unique_genes)) or "None",
        "",
        "Per dataset:",
        "",
    ]

    for dataset_df in dataset_dfs:
        dataset_id: str = str(
            dataset_df.get_column("dataset_id")[0]
        )

        genes: set[str] = genes_per_dataset[dataset_id]

        dataset_unique_genes: set[str] = (
            genes & unique_genes
        )

        dataset_shared_genes: set[str] = (
            genes & shared_genes
        )

        unique_genes_df: pl.DataFrame = (
            dataset_df
            .select(
                "gene_symbol",
                "has_clean_ortholog",
            )
            .drop_nulls("gene_symbol")
            .unique()
        )

        clean_ortholog_df: pl.DataFrame = (
            unique_genes_df
            .filter(
                pl.col("has_clean_ortholog")
                .fill_null(False)
            )
        )

        left_out_df: pl.DataFrame = (
            unique_genes_df
            .filter(
                pl.col("has_clean_ortholog")
                .fill_null(False)
                .not_()
            )
        )

        left_out_genes: list[str] = (
            left_out_df
            .get_column("gene_symbol")
            .cast(pl.String)
            .sort()
            .to_list()
        )

        report_lines.extend(
            [
                dataset_id,
                f"  Total {pattern} genes: {len(genes)}",
                (
                    f"  Shared {pattern} genes: "
                    f"{len(dataset_shared_genes)}"
                ),
                (
                    f"  Dataset-specific {pattern} genes: "
                    f"{len(dataset_unique_genes)}"
                ),
                (
                    f"  {pattern} genes with clean ortholog: "
                    f"{clean_ortholog_df.height}"
                ),
                (
                    f"  {pattern} genes left out from ortholog analysis: "
                    f"{left_out_df.height}"
                ),
                (
                    "  Dataset-specific genes: "
                    f"{', '.join(sorted(dataset_unique_genes))}"
                    if dataset_unique_genes
                    else "  Dataset-specific genes: None"
                ),
                (
                    "  Left out genes: "
                    f"{', '.join(left_out_genes)}"
                    if left_out_genes
                    else "  Left out genes: None"
                ),
                "",
            ]
        )

    report_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    report_file.write_text(
        "\n".join(report_lines),
        encoding="utf-8",
    )

    print(f"Saved report: {report_file}")