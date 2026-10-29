# Define atlas gene sets

Run `run_pipeline.py` with the project analysis Python environment. It runs
all six scripts in dependency order and stops if a step fails. Individual
scripts can also be launched directly from any working directory.

## Order

1. `_01_define_groups/01_define_physiological_groups.py`: group atlas tissues.
2. `_02_define_tsa/01_run_dynamic_step_classification.py`: classify probes.
3. `_02_define_tsa/02_final_gene_cleanup.py`: resolve gene-level TSA calls.
4. `_03_define_non_tsa/01_define_non_tsa.py`: classify non-TSA probes using
   the dynamic-step output.
5. `_03_define_non_tsa/02_final_gene_cleanup.py`: resolve gene-level non-TSA calls.
6. `_04_match_with_datasets/01_match_with_datasets.py`: attach current
   standardized mTEC-high expression statistics to the TSA list.

The non-TSA branch depends on dynamic-step classification, so it must not
run before step 2. Its source folder was previously numbered before TSA
classification. The duplicate-gene helper is now named
`internals/resolve_duplicated_genes.py`.

## Paths

`config/paths.py` defines the input and output locations. Atlas input comes
from `resources/mouse_gene_atlas/standardized/gse10246_atlas.csv`.
Tissue groups and generated files use the existing directory:
`data/mtec_high/derived_data/tsa_expression/samsom_et_al_tsa_analysis/`.

- `tsa_lists/tissue_groups.json`: tissue-group definitions.
- `tsa_lists/raw/`: grouped atlas.
- `tsa_lists/intermediate/`: dynamic-step probe calls.
- `tsa_lists/standardized/`: TSA list and matched expression table.
- `housekeepers_list/intermediate/`: non-TSA probe calls.
- `housekeepers_list/standardized/`: non-TSA gene list.

The existing data directory names are retained. The non-TSA classification
should not be interpreted as independent validation of housekeeping genes.
Expression matching uses the datasets in `config/mtec_high_datasets.toml`
and reads only `data/mtec_high/standardized/`.

Existing classification thresholds and gene-resolution rules are unchanged.
Running the pipeline replaces its generated outputs, but does not modify
the standardized mTEC input files. Sibling tissue-GSEA scripts are outside
this repair.

## Validation

Syntax, local import paths, existing input columns, atlas tissue-group
membership, and runner dispatch were checked. The runner dispatch check
used mocked subprocesses; the full classifications were not rerun.
