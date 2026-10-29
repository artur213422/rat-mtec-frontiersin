from pathlib import Path
import polars as pl

from src.common.paths import PROJECT_ROOT
from src.data_preprocessing.mtec_high.internal.orthology_mapping import OrthologSpecies, find_ortholog_id
from src.resource_processing.michelson_2022_et_al_mimetic_cells.config import p_adj, log2_fold_change, pct_1, pct_2, drop_populations

##########################################################################################
# PREPARE DATA
##########################################################################################

mimetic_cells_file: Path = PROJECT_ROOT / "resources/michelson_2022_et_al_mimetic_cells/source/mimetic_cells.csv"
mimetic_cells_df: pl.DataFrame = pl.read_csv(mimetic_cells_file)
ortholog_map: pl.DataFrame = pl.read_csv(PROJECT_ROOT / "resources/gene_name_conversion/orthologs/HOM_AllOrganism.rpt.csv")

save_file: Path = PROJECT_ROOT / "resources/michelson_2022_et_al_mimetic_cells/standardized/mimetic_cells.gmt"
save_ortholog_file: Path = PROJECT_ROOT / "resources/michelson_2022_et_al_mimetic_cells/standardized/mimetic_cells_ortholog.gmt"

save_summary: Path = PROJECT_ROOT / "resources/michelson_2022_et_al_mimetic_cells/standardized/save_summary.txt"

##########################################################################################
# REGULAR GMT
##########################################################################################

filtered_df: pl.DataFrame = mimetic_cells_df.filter(
    (pl.col("p_adj") < p_adj)
    & (pl.col("log2_fold_change") > log2_fold_change)
    & (pl.col("pct_1") >= pct_1)
    & (pl.col("pct_2") <= pct_2)
)

for pop in drop_populations:
    filtered_df = filtered_df.filter(pl.col("population") != pop)

##########################################################################################
# FIND ORTHOLOG IDS
##########################################################################################

mapping_rows: list[dict[str, str | None]] = []
for gene_symbol in filtered_df.get_column("gene_symbol").drop_nulls().unique().to_list():
    mapping_rows.append({
        "gene_symbol": gene_symbol,
        "ortholog_id": find_ortholog_id(
            gene_symbol,
            OrthologSpecies.MOUSE,
            ortholog_map,
        ),
    })

ortholog_mapping_df = pl.DataFrame(
    mapping_rows,
    schema={
        "gene_symbol": pl.String,
        "ortholog_id": pl.String,
    },
)

filtered_with_orthologs_df = filtered_df.join(
    ortholog_mapping_df,
    on="gene_symbol",
    how="left",
    validate="m:1",
)


##########################################################################################
# SAVE
##########################################################################################

gmt_df: pl.DataFrame = (
    filtered_df
    .group_by("population")
    .agg(
        pl.col("gene_symbol")
        .drop_nulls()
        .unique()
        .sort()
        .alias("genes")
    )
    .sort("population")
)

def write_gmt(gmt_data: pl.DataFrame, output_file: Path) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with output_file.open("w", encoding="utf-8") as file:
        for row in gmt_data.iter_rows(named=True):
            population: str = row["population"]
            genes: list[str] = row["genes"]

            file.write("\t".join([population, "Michelson_2022", *genes]) + "\n")


write_gmt(gmt_df, save_file)

ortholog_ids_df: pl.DataFrame = (
    filtered_with_orthologs_df
    .drop_nulls(subset=["ortholog_id"])
    # A source gene can belong to multiple orthology groups.  Write each group
    # as an individual GMT member; rows without an ortholog were dropped above.
    .with_columns(pl.col("ortholog_id").str.split("///"))
    .explode("ortholog_id")
    .with_columns(pl.col("ortholog_id").str.strip_chars())
    .filter(pl.col("ortholog_id") != "")
)

ortholog_gmt_df: pl.DataFrame = (
    ortholog_ids_df
    .group_by("population")
    .agg(pl.col("ortholog_id").unique().sort().alias("genes"))
    .sort("population")
)

write_gmt(ortholog_gmt_df, save_ortholog_file)


##########################################################################################
# CREATE SUMMARY
##########################################################################################

filtering_summary_df: pl.DataFrame = (
    mimetic_cells_df
    .group_by("population")
    .agg(pl.len().alias("genes_before"))
    .join(
        filtered_with_orthologs_df
        .group_by("population")
        .agg(pl.len().alias("genes_after")),
        on="population",
        how="left",
    )
    .with_columns(
        pl.col("genes_after").fill_null(0),
    )
    .join(
        filtered_with_orthologs_df
        .drop_nulls(subset=["ortholog_id"])
        .group_by("population")
        .agg(pl.len().alias("genes_with_ortholog_id")),
        on="population",
        how="left",
    )
    .join(
        ortholog_ids_df
        .group_by("population")
        .agg(pl.col("ortholog_id").n_unique().alias("unique_ortholog_ids")),
        on="population",
        how="left",
    )
    .with_columns(
        pl.col("genes_with_ortholog_id").fill_null(0),
        pl.col("unique_ortholog_ids").fill_null(0),
    )
    .with_columns(
        (pl.col("genes_after") / pl.col("genes_before") * 100).round(1).alias("retained_pct"),
        (pl.col("genes_with_ortholog_id") / pl.col("genes_after") * 100)
        .fill_nan(0)
        .round(1)
        .alias("ortholog_mapping_pct"),
    )
    .sort("population")
)

# Counts are based on the unique members actually written to each GMT file,
# rather than on the (potentially duplicated) filtered input rows.
gmt_summary_df: pl.DataFrame = (
    gmt_df
    .select(
        "population",
        pl.col("genes").list.len().alias("regular_gmt_gene_count"),
    )
    .join(
        ortholog_gmt_df.select(
            "population",
            pl.col("genes").list.len().alias("ortholog_gmt_id_count"),
        ),
        on="population",
        how="left",
    )
    .with_columns(pl.col("ortholog_gmt_id_count").fill_null(0))
    .sort("population")
)

with pl.Config(tbl_rows=-1, tbl_cols=-1, fmt_str_lengths=100):
    filtering_summary_table: str = str(filtering_summary_df)
    gmt_summary_table: str = str(gmt_summary_df)

summary_text: str = (
    f"Filtering settings\n"
    f"------------------\n"
    f"p_val_adj < {p_adj}\n"
    f"avg_log2FC > {log2_fold_change}\n"
    f"pct.1 >= {pct_1}\n"
    f"pct.2 <= {pct_2}\n\n"
    f"Filtering and ortholog mapping per population\n"
    f"--------------------------------------------\n"
    f"{filtering_summary_table}\n\n"
    f"GMT output per population\n"
    f"-------------------------\n"
    f"regular_gmt_gene_count: unique gene symbols in mimetic_cells.gmt\n"
    f"ortholog_gmt_id_count: unique ortholog IDs in mimetic_cells_ortholog.gmt; "
    f"unmapped genes are omitted and multiple IDs are split.\n"
    f"{gmt_summary_table}\n"
)

save_summary.parent.mkdir(parents=True, exist_ok=True)
save_summary.write_text(summary_text, encoding="utf-8")
