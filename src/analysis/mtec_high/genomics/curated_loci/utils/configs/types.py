from dataclasses import dataclass
from enum import StrEnum
from re import Pattern


class LocusDefinitionType(StrEnum):
    GENE_FAMILY = "gene_family"
    GENE_SET = "gene_set"


@dataclass(frozen=True)
class CuratedLocus:
    locus_name: str
    definition_type: LocusDefinitionType

    gene_patterns: tuple[Pattern[str], ...] = ()
    genes: tuple[str, ...] = ()