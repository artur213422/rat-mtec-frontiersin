# Reference-region expression comparison

Run steps 01, 02 and 03 in order after changing the gene eligibility or reference.
All output paths resolve inside the project containing these scripts.

Step 01 retains the all-gene clustering as a descriptive output, but the
cross-species analysis uses `ortholog_clusters.csv` for BOTH species. Eligibility
requires `has_clean_ortholog=true`, a nonempty ID, and no `///`. Duplicate eligible
IDs within a dataset raise an error. Filtering happens before defining clusters.

The 200 kb rule is an inter-gene gap rule, not a maximum region length. Reference
singletons remain unclustered. Other datasets' singleton genes can map to an
existing reference cluster.

Step 02 assigns mouse genes individually by reference gene symbol. Step 03
assigns rat genes individually by clean ortholog ID. `local_cluster_id` preserves
the original within-dataset region; `matched_cluster_id` identifies the mouse
reference region. One local cluster may contribute to multiple reference regions.
These assignments are not claims of one-to-one conserved physical clusters.

Step 03 outputs:

- `matched_clusters.csv`: shared ortholog pairs only. Rat and reference mouse
  contain the same ortholog set in every included region. Other mouse arrays can
  have missing measurements; use the pairwise gene intersection when comparing.
- `all_cluster_assignments.csv`: all eligible genes, including unmapped genes.
- `excluded_cluster_genes.csv`: eligible genes excluded from the shared output,
  including reference singletons and genes lacking a measured counterpart.
- `cluster_coverage.csv`: measured gene count, local cluster count and local
  singleton count per reference region/dataset. Missing rows mean zero coverage.

The audit tables cover eligible clean orthologs; non-clean genes remain available
in the all-gene source table. Mouse matching has its own summary and does not
overwrite the local-clustering summary.

Regression checks: `python -m unittest src.analysis.mtec_high.simple_clustering.test_matching`
