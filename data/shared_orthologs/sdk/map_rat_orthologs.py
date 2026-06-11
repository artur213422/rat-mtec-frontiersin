import time
from pathlib import Path
import polars as pl

from data.shared_orthologs.sdk.helpers.helpers import (
    OrthologResult,
    Species,
    find_ortholog,
    Status,
)

mtec_high_file = Path(r"G:\10. ÜLIÕPILASED\Artur\9. Review\analysis\data\rat\cleaned_data\cleaned_mtec_high.csv")
mtec_low_file = Path(r"G:\10. ÜLIÕPILASED\Artur\9. Review\analysis\data\rat\cleaned_data\cleaned_mtec_low.csv")
ortholog_map_file = Path(r"G:\10. ÜLIÕPILASED\Artur\9. Review\analysis\data\shared_orthologs\HOM_AllOrganism.rpt.csv")

output = Path(r"G:\10. ÜLIÕPILASED\Artur\9. Review\analysis\data\shared_orthologs\mapped\rat")
output.mkdir(parents=True, exist_ok=True)

analyse_file: Path = mtec_high_file

mtec_df: pl.DataFrame = pl.read_csv(analyse_file)
ortholog_map: pl.DataFrame = pl.read_csv(ortholog_map_file)

ok_count = 0
no_orthology_count = 0
no_orthology_id_count = 0
rat_polymorphism_count = 0
mouse_polymorphism_count = 0

mapped_rows = []

start_time = time.perf_counter()

for row in mtec_df.iter_rows(named=True):
    gene: str = row["gene_symbol"]
    result: OrthologResult = find_ortholog(gene, Species.RAT, ortholog_map)

    match result.status:
        case Status.OK:
            ok_count += 1

            if len(result.rat_symbols) > 1:
                rat_polymorphism_count += 1

            if len(result.mouse_symbols) > 1:
                mouse_polymorphism_count += 1

        case Status.NO_ORTHOLOGY_FOUND:
            no_orthology_count += 1

        case Status.ORTHOLOGY_ID_NOT_FOUND:
            no_orthology_id_count += 1

    mapped_rows.append({
        **row,

        # new fields
        "ortholog_status": result.status.value,
        "ortholog_id": result.ortholog_id,
        "mouse_ortholog_name": ";".join(result.mouse_symbols),

        # optional but useful audit fields
        "rat_ortholog_names": ";".join(result.rat_symbols),
        "n_mouse_orthologs": len(result.mouse_symbols),
        "n_rat_orthologs": len(result.rat_symbols),
    })

elapsed = time.perf_counter() - start_time

mapped_df = pl.DataFrame(mapped_rows)

clean_name = analyse_file.stem.removeprefix("cleaned_")
output_file = output / f"{clean_name}_mapped_mouse_orthologs.csv"

mapped_df.write_csv(output_file)

print("---- Summary ----")
print(mtec_df.shape)
print(mtec_df.height)
print(f"OK: {ok_count}")
print(f"No orthology found: {no_orthology_count}")
print(f"Orthology ID not found: {no_orthology_id_count}")

print(f"Polymorphic mouse genes: {mouse_polymorphism_count}")
print(f"Polymorphic rat genes: {rat_polymorphism_count}")
print(f"Total: {ok_count + no_orthology_count + no_orthology_id_count}")
print(f"Time elapsed: {elapsed:.2f} seconds")

print(f"Saved mapped file to: {output_file}")