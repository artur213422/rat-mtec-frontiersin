from pathlib import Path

from src.common.paths import PROJECT_ROOT
from src.resource_processing.standard_gsea_pathways.internal.enums import Namespace


CORE_PATHWAYS: dict[Namespace, list[str]] = {
    Namespace.BP: [
        "GO:0022613 | ribonucleoprotein complex biogenesis",
        "GO:0042254 | ribosome biogenesis",
        "GO:0016072 | rRNA metabolic process",
        "GO:0006364 | rRNA processing",
        "GO:0006260 | DNA replication",
        "GO:0006275 | regulation of DNA replication",
        "GO:0006310 | DNA recombination",
        "GO:0006302 | double-strand break repair",
        "GO:0006405 | RNA export from nucleus",
        "GO:0000245 | spliceosomal complex assembly",
    ],
    Namespace.CP: [
        "GO:0030684 | preribosome",
        "GO:0030686 | 90S preribosome",
        "GO:0030687 | preribosome, large subunit precursor",
        "GO:0032040 | small-subunit processome",
    ],
    Namespace.MF: [],
}


GSEA_FOLDER: Path = PROJECT_ROOT / "data/mtec_high/derived_data/standard_gsea/both_sided_gsea"

GSEA_ANALYSIS: dict[Namespace, Path] = {
    Namespace.BP: GSEA_FOLDER / "go_biological_process_gsea.csv",
    Namespace.CP: GSEA_FOLDER / "go_cellular_component_gsea.csv",
    Namespace.MF: GSEA_FOLDER / "go_molecular_function_gsea.csv",
}

MIN_DATASET_SUPPORT: int = 2


