import polars as pl


def resolve_duplicated_genes(atlas_df: pl.DataFrame, maximum_selected_groups: int = 5) -> pl.DataFrame:
    required_columns: set[str] = {
        "ID_REF",
        "ortholog_id",
        "gene_symbol",
        "is_tissue_restricted",
        "selected_groups",
    }

    missing_columns: set[str] = required_columns - set(atlas_df.columns)

    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    if maximum_selected_groups < 1:
        raise ValueError("maximum_selected_groups must be at least 1")

    if atlas_df["ortholog_id"].null_count() > 0:
        raise ValueError("Input contains null ortholog IDs")

    resolved_rows: list[dict] = []

    gene_tables: list[pl.DataFrame] = atlas_df.partition_by(
        "ortholog_id",
        maintain_order=True,
    )

    ######################################################################################
    # REPORT COUNTERS
    ######################################################################################

    total_gene_count: int = len(gene_tables)

    no_restricted_probe_count: int = 0
    unanimously_restricted_count: int = 0

    duplicated_gene_count: int = 0
    differing_tissue_call_count: int = 0

    any_single_call_became_multi_count: int = 0
    all_single_calls_became_multi_count: int = 0

    dropped_too_many_groups_count: int = 0

    ######################################################################################
    # RESOLVE GENES
    ######################################################################################

    supported_gene_count = 0

    for gene_df in gene_tables:
        ortholog_id = gene_df.item(0, "ortholog_id")
        probe_count = gene_df.height

        if probe_count > 1:
            duplicated_gene_count += 1

        restriction_calls: list[bool] = (
            gene_df["is_tissue_restricted"]
            .fill_null(False)
            .to_list()
        )

        # Require at least one tissue-restricted probe set.
        if not any(restriction_calls):
            no_restricted_probe_count += 1
            continue

        supported_gene_count += 1

        gene_symbols: list[str] = sorted({
            symbol.strip()
            for symbol in gene_df["gene_symbol"].to_list()
            if symbol is not None and symbol.strip()
        })

        if not gene_symbols:
            raise ValueError(
                f"No valid gene symbol for ortholog {ortholog_id}"
            )

        probe_ids: list[str] = sorted({
            probe_id
            for probe_id in gene_df["ID_REF"].to_list()
            if probe_id is not None
        })

        normalized_probe_calls: list[tuple[str, ...]] = []

        for is_restricted, selected_groups in zip(restriction_calls, gene_df["selected_groups"].to_list()):
            if not is_restricted or selected_groups is None:
                continue

            tissues: tuple[str, ...] = tuple(sorted({
                tissue.strip()
                for tissue in selected_groups.split(";")
                if tissue.strip()
            }))

            if tissues:
                normalized_probe_calls.append(tissues)

        if not normalized_probe_calls:
            raise ValueError(f"Unanimously restricted gene has no selected tissues: ortholog {ortholog_id}")

        unique_probe_calls: set[tuple[str, ...]] = set(normalized_probe_calls)

        if len(unique_probe_calls) > 1:
            differing_tissue_call_count += 1

        # Union tissue assignments from all probes.
        combined_tissues: set[str] = {
            tissue
            for probe_call in normalized_probe_calls
            for tissue in probe_call
        }

        final_primary_j: int = len(combined_tissues)

        original_primary_j_values: list[int] = [
            len(probe_call)
            for probe_call in normalized_probe_calls
        ]

        # At least one probe was single-tissue, but union became multi-tissue.
        if (
            final_primary_j > 1
            and any(
                primary_j == 1
                for primary_j in original_primary_j_values
            )
        ):
            any_single_call_became_multi_count += 1

        # Every probe was individually single-tissue, but probes selected
        # different tissues and the union became multi-tissue.
        if (
            final_primary_j > 1
            and all(
                primary_j == 1
                for primary_j in original_primary_j_values
            )
        ):
            all_single_calls_became_multi_count += 1

        # A union broader than the original dynamic-step limit is invalid.
        if final_primary_j > maximum_selected_groups:
            dropped_too_many_groups_count += 1
            continue

        resolved_rows.append(
            {
                "ortholog_id": ortholog_id,
                "gene_symbol": gene_symbols[0],
                "is_tissue_restricted": True,
                "primary_j": final_primary_j,
                "selected_groups": ";".join(
                    sorted(combined_tissues)
                ),
                "probe_count": len(probe_ids),
                "probe_calls_agree": (
                    len(unique_probe_calls) == 1
                ),
                "gene_symbol_conflict": (
                    len(gene_symbols) > 1
                ),
                "probe_ids": ";".join(probe_ids),
                "probe_group_calls": " | ".join(
                    ";".join(probe_call)
                    for probe_call in sorted(
                        unique_probe_calls
                    )
                ),
            }
        )

    ######################################################################################
    # BUILD OUTPUT
    ######################################################################################

    resolved_df: pl.DataFrame = pl.DataFrame(
        resolved_rows,
        schema={
            "ortholog_id": pl.Int64,
            "gene_symbol": pl.String,
            "is_tissue_restricted": pl.Boolean,
            "primary_j": pl.Int64,
            "selected_groups": pl.String,
            "probe_count": pl.Int64,
            "probe_calls_agree": pl.Boolean,
            "gene_symbol_conflict": pl.Boolean,
            "probe_ids": pl.String,
            "probe_group_calls": pl.String,
        },
        strict=False,
    )

    final_single_tissue_count: int = resolved_df.filter(pl.col("primary_j") == 1).height

    final_multi_tissue_count: int = resolved_df.filter(pl.col("primary_j") > 1).height

    ######################################################################################
    # PRINT REPORT
    ######################################################################################

    print()
    print("PROBE-TO-GENE TSA RESOLUTION")
    print("=" * 50)
    print(f"Total ortholog genes examined: {total_gene_count}")
    print(f"Genes represented by multiple probes: {duplicated_gene_count}")
    print()
    print(f"NRE genes dropped: {no_restricted_probe_count}")
    print(f"Unanimously tissue-restricted genes: {unanimously_restricted_count}")
    print()
    print(f"Genes with differing restricted-probe tissue calls: {differing_tissue_call_count}")
    print(f"Genes where at least one primary_j=1 probe became multi-tissue after union: {any_single_call_became_multi_count}")
    print(f"Genes where all probes had primary_j=1 but their union became multi-tissue: {all_single_calls_became_multi_count}")
    print(f"Genes dropped because union exceeded {maximum_selected_groups} groups: {dropped_too_many_groups_count}")
    print()
    print(f"Final retained genes: {resolved_df.height}")
    print(f"Final single-tissue genes: {final_single_tissue_count}")
    print(f"Final multi-tissue genes: {final_multi_tissue_count}")
    print("=" * 50)

    return resolved_df