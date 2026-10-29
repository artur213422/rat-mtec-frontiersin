"""Input-validated, atomic checkpoints for completed GSEA calls."""
from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path

from src.common.feratures.gsea.gsea import run_gsea
from src.common.feratures.gsea.models import GseaInput, GseaResult


def run_gsea_checkpointed(
    input: GseaInput, *, checkpoint_file: Path, dataset_id: str, species: str,
) -> list[GseaResult]:
    source_folder = Path(__file__).resolve().parent
    payload = {
        "version": 1,
        "dataset_id": dataset_id,
        "species": species,
        "ranked_genes": input.rnk.to_dicts(),
        "gene_sets": input.gene_sets_dict,
        "min_size": input.min_size,
        "max_size": input.max_size,
        "permutation_num": input.permutation_num,
        "mode": input.gsea_mode.value,
        "implementation": {
            name: sha256((source_folder / name).read_bytes()).hexdigest()
            for name in ("gsea.R", "gsea.py", "models.py", "checkpoint.py")
        },
    }
    fingerprint = sha256(json.dumps(payload, sort_keys=True, allow_nan=False).encode()).hexdigest()
    if checkpoint_file.is_file():
        try:
            saved = json.loads(checkpoint_file.read_text(encoding="utf-8"))
            if saved["fingerprint"] == fingerprint:
                results = [GseaResult(**row) for row in saved["results"]]
                print(f"Resuming completed GSEA: {checkpoint_file.name}")
                return results
        except (ValueError, KeyError, TypeError):
            print(f"Ignoring invalid checkpoint: {checkpoint_file}")

    results = run_gsea(input)
    checkpoint_file.parent.mkdir(parents=True, exist_ok=True)
    temporary_file = checkpoint_file.with_suffix(".tmp")
    temporary_file.write_text(json.dumps({
        "fingerprint": fingerprint,
        "dataset_id": dataset_id,
        "species": species,
        "mode": input.gsea_mode.value,
        "results": [asdict(result) for result in results],
    }, allow_nan=False), encoding="utf-8")
    temporary_file.replace(checkpoint_file)
    return results