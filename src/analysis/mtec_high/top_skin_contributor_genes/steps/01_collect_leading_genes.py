from pathlib import Path

import polars as pl

from src.analysis.mtec_high.genomics.window_sliding.run_pipeline import PROJECT_ROOT
from src.analysis.mtec_high.top_skin_contributor_genes.helpers.collect_candidates import collect_michelson_gsea, \
    collect_hollander_gsea, collect_tsa_gsea, collect_go_gsea

# ------------------------------------------------------------------------------------------------
# Prepare Data
# ------------------------------------------------------------------------------------------------

output_file: Path = PROJECT_ROOT / "data/mtec_high/derived_data/top_skin_contributor_genes/source/01_skin_leading_genes.csv"

# ------------------------------------------------------------------------------------------------
# Collect Candidates
# ------------------------------------------------------------------------------------------------

candidates: pl.DataFrame = pl.DataFrame()

candidates = collect_go_gsea(candidates)
candidates = collect_tsa_gsea(candidates)
candidates = collect_hollander_gsea(candidates)
candidates = collect_michelson_gsea(candidates)

# ------------------------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------------------------

candidates.write_csv(output_file)


