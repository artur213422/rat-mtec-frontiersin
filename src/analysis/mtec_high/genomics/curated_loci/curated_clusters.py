import re
from re import Pattern

from src.analysis.mtec_high.genomics.curated_loci.utils.configs.types import LocusDefinitionType, CuratedLocus

import re

from src.analysis.mtec_high.genomics.curated_loci.utils.configs.types import LocusDefinitionType, CuratedLocus

GENE_SETS: dict[str, tuple[str, ...]] = {
    "mhc": (
        "H2-K1",
        "H2-D1",
        "H2-Aa",
        "H2-Ab1",
        "Tap1",
        "Tap2",
    ),
    "insulin": (
        "Ins2",
    ),
    "keratin_type_i": (
        "Krt14",
        "Krt16",
        "Krt17",
    ),
    "keratin_type_ii": (
        "Krt5",
        "Krt6a",
        "Krt6b",
    ),
    "salivary_protein": (
        "Spt1",
    ),
}


GENE_FAMILIES: dict[str, CuratedLocus] = {
    # Family membership is a candidate-selection rule, not evidence that all
    # members share one genomic locus. Inspect coordinates in each species.
    "casein": CuratedLocus(
        locus_name="casein",
        definition_type=LocusDefinitionType.GENE_FAMILY,
        gene_patterns=(
            re.compile(r"^Csn(?:1s[12][a-z]?|2|3)$", re.IGNORECASE),
        ),
    ),
    "major_urinary_protein": CuratedLocus(
        locus_name="major_urinary_protein",
        definition_type=LocusDefinitionType.GENE_FAMILY,
        gene_patterns=(
            re.compile(r"^Mup\d+[a-z]?$", re.IGNORECASE),
        ),
    ),
    "mucin": CuratedLocus(
        locus_name="mucin",
        definition_type=LocusDefinitionType.GENE_FAMILY,
        gene_patterns=(
            re.compile(r"^Muc\d+[a-z]*$", re.IGNORECASE),
        ),
    ),
    "kallikrein": CuratedLocus(
        locus_name="kallikrein",
        definition_type=LocusDefinitionType.GENE_FAMILY,
        gene_patterns=(
            re.compile(r"^Klk(?:[1-9]|1[0-5])$", re.IGNORECASE),
            re.compile(r"^Klk1b\d+[a-z]?$",re.IGNORECASE),
            re.compile(r"^Klkb1$", re.IGNORECASE),
            re.compile(r"^Klk1c\d+[a-z]?$", re.IGNORECASE),
            re.compile(r"^Klk5l$", re.IGNORECASE),
        ),
    ),
}
