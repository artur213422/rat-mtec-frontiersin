"""Explicit inputs and outputs for the atlas gene-set pipeline."""


from src.common.paths import PROJECT_ROOT

DATA_FOLDER = PROJECT_ROOT / "data/mtec_high/derived_data/tsa_expression/samsom_et_al_tsa_analysis"
STANDARDIZED_FOLDER = PROJECT_ROOT / "data/mtec_high/standardized"
ATLAS_FILE = PROJECT_ROOT / "resources/mouse_gene_atlas/standardized/gse10246_atlas.csv"
TISSUE_GROUP_FILE = DATA_FOLDER / "tsa_lists/tissue_groups.json"
GROUPED_ATLAS_FILE = DATA_FOLDER / "tsa_lists/raw/mouse_gene_atlas_tissue_groups.csv"
DYNAMIC_STEP_FILE = DATA_FOLDER / "tsa_lists/intermediate/gse10246_dynamic_step_probe_calls.csv"
TSA_FILE = DATA_FOLDER / "tsa_lists/standardized/mouse_tissue_atlas.csv"
NON_TSA_PROBE_FILE = DATA_FOLDER / "housekeepers_list/intermediate/non_tsa_probe_calls.csv"
NON_TSA_FILE = DATA_FOLDER / "housekeepers_list/standardized/non_tsa_list.csv"
MATCHED_TSA_FILE = DATA_FOLDER / "tsa_lists/standardized/tsa_lists_with_rodent_expression.csv"
