# ISG pathway analysis

Outputs are saved under `data/mtec_high/derived_data/isg_pathways/`.
Inputs come from `resources/interferon/`, the shared ortholog resource,
and `data/mtec_high/standardized/`.

Run each branch in numerical order. Step 01 creates human-gene pathway
statistics; step 02 maps those genes to mouse orthologs in the same CSV.
When rerunning ortholog mapping on older results that lack the original
`gene` column, rerun step 01 first.

The lung IFN-versus-all step 03 optionally exports clean-ortholog gene lists
and a CSV summary under `isg_pathway_genes/picked_for_gsea/`. Significant
upregulated genes were already selected in step 01. Step 04 continues to
use the mapped table, including fold changes needed for its IFN-beta
500-gene limit; step 03 is not a prerequisite for GSEA.

Existing ranking calculations and branch-specific IFN gene limits are
preserved. Human-TEC GSEA uses standardized ranking statistics, while the
lung branches retain their original ranking calculations.

Use the project analysis environment, including Polars, GSEA dependencies,
and the configured R/fgsea runtime. The scripts support direct execution
from any working directory.
