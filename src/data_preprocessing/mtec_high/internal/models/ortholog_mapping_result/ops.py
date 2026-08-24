from dataclasses import fields
from pathlib import Path

from src.data_preprocessing.mtec_high.internal.models.ortholog_mapping_result.model import OrthologMappingResult


def save_ortholog_mapping_result(result: OrthologMappingResult, output_path: Path, filename: str) -> None:
    output_path.mkdir(parents=True, exist_ok=True)

    output_file: Path = output_path / f"{filename}_ortholog_mapping_result.txt"

    with output_file.open("w", encoding="utf-8") as file:
        for field in fields(result):
            value = getattr(result, field.name)
            file.write(f"{field.name}: {value}\n")