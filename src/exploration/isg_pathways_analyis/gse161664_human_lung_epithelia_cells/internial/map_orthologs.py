import polars as pl

from src.data_preprocessing.internal.orthology_mapping import OrthologRequest, find_ortholog, OrthologSpecies, Status


def map_orthologs(interferon_isg: pl.DataFrame, ortholog_map: pl.DataFrame) -> pl.DataFrame:
    ortholog_ids: list[str | None] = []
    mouse_symbols: list[str | None] = []
    has_clean_orthologs: list[bool] = []

    for row in interferon_isg.iter_rows(named=True):
        gene: str = row["gene"]

        result: OrthologRequest = find_ortholog(
            gene,
            OrthologSpecies.HUMAN,
            ortholog_map,
        )

        has_clean_ortholog: bool = True

        match result.status:
            case Status.OK:
                if len(result.rat_symbols) > 1:
                    has_clean_ortholog = False
                if len(result.mouse_symbols) > 1:
                    has_clean_ortholog = False

            case Status.ORTHOLOGY_ID_NOT_FOUND:
                has_clean_ortholog = False

        gene_symbol: str = ", ".join(result.mouse_symbols)

        ortholog_ids.append(result.ortholog_id)
        has_clean_orthologs.append(has_clean_ortholog)
        mouse_symbols.append(gene_symbol)

    interferon_isg = interferon_isg.with_columns(
        pl.Series("ortholog_id", ortholog_ids),
        pl.Series("gene_symbol", mouse_symbols),
        pl.Series("has_clean_ortholog", has_clean_orthologs),
    )


    return interferon_isg
