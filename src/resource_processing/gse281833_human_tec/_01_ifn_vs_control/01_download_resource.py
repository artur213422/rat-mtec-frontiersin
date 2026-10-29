from pathlib import Path

from src.common.file_io import (
    delete_file,
    download_from_web,
    unpack_file,
)
from src.common.paths import resolve_dir_from_root

GSE281833_URL: str = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281833/suppl/GSE281833_count_matrix.tsv.gz"


output_dir: Path = resolve_dir_from_root("resources/interferon/gse28183_tec/source/")

def download_resource() -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)

    archive_path: Path = output_dir / "GSE281833_count_matrix.tsv.gz"
    output_path: Path = output_dir / "GSE281833_count_matrix.tsv"

    download_from_web(GSE281833_URL, archive_path)

    unpack_file(archive_path, output_dir)

    delete_file(archive_path)

    return output_path


count_matrix_path: Path = download_resource()
print(count_matrix_path)
