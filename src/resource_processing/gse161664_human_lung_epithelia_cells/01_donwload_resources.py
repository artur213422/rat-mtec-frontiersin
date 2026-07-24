from pathlib import Path

from src.common.file_io import (
    delete_file,
    download_from_web,
    unpack_file,
)
from src.common.paths import resolve_dir_from_root


GSE161664_RAW_URL: str = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161664/suppl/GSE161664_RAW.tar"
)

output_dir: Path = resolve_dir_from_root("resources/interferon/gse161664")

def download_gse161664_resources() -> None:


    archive_path: Path = output_dir / "GSE161664_RAW.tar"

    download_from_web(
        GSE161664_RAW_URL,
        archive_path,
    )

    unpack_file(
        archive_path,
        output_dir,
    )

    delete_file(archive_path)

download_gse161664_resources()
paths: list[Path] = list(output_dir.glob("*"))

for file in paths:
    unpack_file(file, output_dir)
    delete_file(file)
