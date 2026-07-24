import polars as pl
import re

from src.common.models.dataset_metadata.model import DatasetMetadata
from src.common.types.enums import Genotypes
from src.common.utils.col_names import GENE_SYMBOL_COLUMN, ORTHOLOG_ID_COLUMN


def _get_genotype_columns(df: pl.DataFrame, gt: Genotypes) -> list[str]:
    pattern: re.Pattern[str] = re.compile(rf"^{re.escape(gt.value)}\d+$")

    return [
        column
        for column in df.columns
        if pattern.fullmatch(column) is not None
    ]


def _get_clean_ortholog_filter(df: pl.DataFrame) -> pl.Expr:
    if "has_clean_ortholog" in df.columns:
        column: str = "has_clean_ortholog"
    elif "has_ortholog" in df.columns:
        column = "has_ortholog"
    else:
        raise ValueError("Missing clean ortholog column. Expected 'has_clean_ortholog' or 'has_ortholog'.")

    if df.schema[column] == pl.Boolean:
        return pl.col(column)

    return (
        pl.col(column)
        .cast(pl.Utf8)
        .str.to_lowercase()
        .is_in(["true", "1", "yes"])
    )


def _read_dataset_df(metadata: DatasetMetadata,use_clean_orthologs: bool) -> pl.DataFrame:
    df: pl.DataFrame = metadata.get_standardized_df()

    if use_clean_orthologs:
        return df.filter(_get_clean_ortholog_filter(df))

    return df




def _build_dataset_frame(
    metadata: DatasetMetadata,
    join_columns: tuple[str, ...],
    ctrl_start_index: int,
    ko_start_index: int,
    use_clean_orthologs: bool,
    keep_gene_symbol_audit: bool,
    keep_ortholog_id_audit: bool,
) -> tuple[
    pl.DataFrame,
    int,
    int,
    list[str],
    list[str],
    list[str],
    list[str],
]:
    df: pl.DataFrame = _read_dataset_df(
        metadata=metadata,
        use_clean_orthologs=use_clean_orthologs,
    )

    missing_id_columns: list[str] = [
        column
        for column in join_columns
        if column not in df.columns
    ]

    if missing_id_columns:
        raise ValueError(
            f"{metadata.dataset_id}: missing columns: {missing_id_columns}"
        )

    ctrl_cols: list[str] = _get_genotype_columns(df, Genotypes.CTRL)
    ko_cols: list[str] = _get_genotype_columns(df, Genotypes.KO)

    if not ctrl_cols:
        raise ValueError(f"{metadata.dataset_id}: no CTRL columns found")

    if not ko_cols:
        raise ValueError(f"{metadata.dataset_id}: no KO columns found")

    rename_map: dict[str, str] = {}

    new_ctrl_cols: list[str] = []
    for old_col in ctrl_cols:
        new_col: str = f"ctrl{ctrl_start_index}"
        rename_map[old_col] = new_col
        new_ctrl_cols.append(new_col)
        ctrl_start_index += 1

    new_ko_cols: list[str] = []
    for old_col in ko_cols:
        new_col: str = f"ko{ko_start_index}"
        rename_map[old_col] = new_col
        new_ko_cols.append(new_col)
        ko_start_index += 1

    selected_columns: list[str] = [
        *join_columns,
        *ctrl_cols,
        *ko_cols,
    ]

    audit_symbol_cols: list[str] = []
    audit_ortholog_cols: list[str] = []

    if (
        keep_gene_symbol_audit
        and GENE_SYMBOL_COLUMN in df.columns
        and GENE_SYMBOL_COLUMN not in join_columns
    ):
        audit_gene_symbol_col: str = f"{metadata.dataset_id}_gene_symbol"
        selected_columns.append(GENE_SYMBOL_COLUMN)
        rename_map[GENE_SYMBOL_COLUMN] = audit_gene_symbol_col
        audit_symbol_cols.append(audit_gene_symbol_col)

    if (
        keep_ortholog_id_audit
        and ORTHOLOG_ID_COLUMN in df.columns
        and ORTHOLOG_ID_COLUMN not in join_columns
    ):
        audit_ortholog_col: str = f"{metadata.dataset_id}_ortholog_id"
        selected_columns.append(ORTHOLOG_ID_COLUMN)
        rename_map[ORTHOLOG_ID_COLUMN] = audit_ortholog_col
        audit_ortholog_cols.append(audit_ortholog_col)

    selected_df: pl.DataFrame = (
        df
        .select(selected_columns)
        .rename(rename_map)
    )

    return (
        selected_df,
        ctrl_start_index,
        ko_start_index,
        new_ctrl_cols,
        new_ko_cols,
        audit_symbol_cols,
        audit_ortholog_cols,
    )




def _combine_datasets_with_nulls(
    datasets: list[DatasetMetadata],
    join_columns: tuple[str, ...],
    use_clean_orthologs: bool,
    keep_gene_symbol_audit: bool,
    keep_ortholog_id_audit: bool,
) -> pl.DataFrame:
    if not datasets:
        return pl.DataFrame(
            schema={
                column: pl.Utf8
                for column in join_columns
            }
        )

    frames: list[pl.DataFrame] = []

    ctrl_index: int = 1
    ko_index: int = 1

    all_ctrl_cols: list[str] = []
    all_ko_cols: list[str] = []
    all_audit_symbol_cols: list[str] = []
    all_audit_ortholog_cols: list[str] = []

    for metadata in datasets:
        (
            frame,
            ctrl_index,
            ko_index,
            ctrl_cols,
            ko_cols,
            audit_symbol_cols,
            audit_ortholog_cols,
        ) = _build_dataset_frame(
            metadata=metadata,
            join_columns=join_columns,
            ctrl_start_index=ctrl_index,
            ko_start_index=ko_index,
            use_clean_orthologs=use_clean_orthologs,
            keep_gene_symbol_audit=keep_gene_symbol_audit,
            keep_ortholog_id_audit=keep_ortholog_id_audit,
        )

        frames.append(frame)
        all_ctrl_cols.extend(ctrl_cols)
        all_ko_cols.extend(ko_cols)
        all_audit_symbol_cols.extend(audit_symbol_cols)
        all_audit_ortholog_cols.extend(audit_ortholog_cols)

    combined: pl.DataFrame = frames[0]

    for frame in frames[1:]:
        combined = combined.join(
            frame,
            on=list(join_columns),
            how="full",
            coalesce=True,
        )

    if all_audit_symbol_cols:
        combined = combined.with_columns(
            pl.coalesce(
                [pl.col(column).cast(pl.Utf8) for column in all_audit_symbol_cols]
            ).alias(GENE_SYMBOL_COLUMN)
        )

    if all_audit_ortholog_cols:
        combined = combined.with_columns(
            pl.coalesce(
                [pl.col(column).cast(pl.Utf8) for column in all_audit_ortholog_cols]
            ).alias(ORTHOLOG_ID_COLUMN)
        )

    #print_null_gene_summary(
    #    data=combined,
   #     id_columns=join_columns,
    #)

    final_columns: list[str] = [
        *join_columns,
    ]

    if (
        GENE_SYMBOL_COLUMN in combined.columns
        and GENE_SYMBOL_COLUMN not in final_columns
    ):
        final_columns.append(GENE_SYMBOL_COLUMN)

    if (
        ORTHOLOG_ID_COLUMN in combined.columns
        and ORTHOLOG_ID_COLUMN not in final_columns
    ):
        final_columns.append(ORTHOLOG_ID_COLUMN)

    final_columns.extend(all_ctrl_cols)
    final_columns.extend(all_ko_cols)

    return combined.select(final_columns)
