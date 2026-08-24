from src.common.paths import find_file_in_project, resolve_dir_from_root
from pathlib import Path

import polars as pl

from src.exploration.isg_pathways_analyis.gse161664_human_lung_epithelia_cells.internial.pick_isg import is_candidate

##########################################################################################
#                                       PREPARE DATA
##########################################################################################


output_file: Path = resolve_dir_from_root("data/analysis/isg_pathways_analysis/gse161664_human_lung_epithelia_cells/isg_patterns/picked_for_gsea") / "gse161664_iss_picked_for_gsea.csv"
summary_file: Path = resolve_dir_from_root("data/analysis/isg_pathways_analysis/gse161664_human_lung_epithelia_cells/isg_patterns/picked_for_gsea") / "summary_file.csv"

source_folder: Path = resolve_dir_from_root("data/analysis/isg_pathways_analysis/gse161664_human_lung_epithelia_cells/isg_patterns/raw")
files: list[Path] = list(source_folder.glob("*.csv"))

##########################################################################################
#                                       PROCESS DATA
##########################################################################################

isg_dfs: list[pl.DataFrame] = []
summary_rows: list[dict[str, str | int]] = []


for file in files:
    interferon_isg: pl.DataFrame = pl.read_csv(file)

    mouse_symbols: list[str] = []
    pathways: list[str] = []

    picked_count: int = 0
    skipped_count: int = 0

    skipped_no_ortholog: int = 0
    skipped_not_candidate: int = 0
    skipped_invalid_symbol: int = 0

    for row in interferon_isg.iter_rows(named=True):
        if row["has_clean_ortholog"] is False:
            skipped_count += 1
            skipped_no_ortholog += 1
            continue
        if row["gene_symbol"] is None:
            continue
        if not is_candidate(row):
            skipped_count += 1
            skipped_not_candidate += 1
            continue

        mouse_symbol: object = row["gene_symbol"]

        if not isinstance(mouse_symbol, str):
            skipped_count += 1
            skipped_invalid_symbol += 1
            continue

        mouse_symbols.append(mouse_symbol)
        pathways.append(file.stem)

        picked_count += 1

    picked_df: pl.DataFrame = pl.DataFrame(
        {
            "gene_symbol": mouse_symbols,
            "pathway": pathways,
        }
    )

    isg_dfs.append(picked_df)

    summary_rows.append(
        {
            "pathway": file.stem,
            "total_genes": interferon_isg.height,
            "picked_genes": picked_count,
            "skipped_genes": skipped_count,
            "skipped_no_ortholog": skipped_no_ortholog,
            "skipped_not_candidate": skipped_not_candidate,
            "skipped_invalid_symbol": skipped_invalid_symbol,
        }
    )


isg_df: pl.DataFrame = pl.concat(
    isg_dfs,
    how="vertical",
)

summary_df: pl.DataFrame = pl.DataFrame(summary_rows)


##########################################################################################
#                                           SAVE
##########################################################################################

isg_df.write_csv(output_file)
summary_df.write_csv(summary_file)

print(isg_df)
print(summary_df)