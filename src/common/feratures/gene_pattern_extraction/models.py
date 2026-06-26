from dataclasses import dataclass




@dataclass(slots=True)
class GeneReferencePerDataset:
    dataset_id: str
    gen_refs: list[GeneReference]


@dataclass(slots=True)
class SharedGeneReferences:
    missing_refs: list[GeneReference]
    shared_refs: list[GeneReference]