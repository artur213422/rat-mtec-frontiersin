from dataclasses import fields
from pathlib import Path

from src.data_preprocessing.internal.models.preprocessing_result.model import PreprocessingResult


def save_preprocessing_result(result: PreprocessingResult, output_path: Path, filename: str) -> None:
    output_path.mkdir(parents=True, exist_ok=True)

    output_file: Path = output_path / f"{filename}_preprocessing_result.txt"

    with output_file.open("w", encoding="utf-8") as file:
        for field in fields(result):
            value = getattr(result, field.name)
            file.write(f"{field.name}: {value}\n")