from pathlib import Path
from pydoc import resolve

import polars as pl

from src.common.paths import find_file_in_project, resolve_dir_from_root
from src.data_preprocessing.mtec_high.internal.orthology_mapping import find_ortholog, OrthologRequest, OrthologSpecies, \
    Status

##########################################################################################
# PREPARE DATA
##########################################################################################

aps1_antigenes_file: Path = find_file_in_project("aps1_aire_antigens_source.csv")

aps1_antigenes_df: pl.DataFrame = pl.read_csv(aps1_antigenes_file)

ortholog_map_file: Path = find_file_in_project("HOM_AllOrganism.rpt.csv")
ortholog_map: pl.DataFrame = pl.read_csv(ortholog_map_file)

out_file: Path = resolve_dir_from_root("data/mtec_high/derived_data/tsa_expression/aps1_known_antigens/source/intermediate") / "aps1_antigens_with_orthologs.csv"

##########################################################################################
# FIND ORTHOLOG
##########################################################################################


mapped_rows: list[dict[str, object]] = []

ortholog_ids: list[str | None] = []
mouse_symbols: list[str | None] = []
has_clean_orthologs: list[bool] = []

for row in aps1_antigenes_df.iter_rows(named=True):
    gene: str = row["human_gene_symbol"]

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

aps1_antigenes_df = aps1_antigenes_df.with_columns(
    pl.Series("ortholog_id", ortholog_ids),
    pl.Series("gene_symbol", mouse_symbols),
    pl.Series("has_clean_ortholog", has_clean_orthologs),
)




##########################################################################################
# SAVE
##########################################################################################
aps1_antigenes_df.write_csv(out_file)

