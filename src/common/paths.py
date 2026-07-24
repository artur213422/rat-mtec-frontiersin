from pathlib import Path

PROJECT_ROOT: Path = Path(__file__).resolve().parents[2]

def try_find_file_in_project(file_name: str, root: Path = PROJECT_ROOT) -> Path | None:
    matches: list[Path] = [
        path.resolve()
        for path in root.rglob(file_name)
        if path.is_file()
    ]

    if not matches:
        return None

    if len(matches) > 1:
        raise ValueError(f"Found multiple files named {file_name!r}: {[str(path) for path in matches]}")

    return matches[0]


def find_file_in_project(file_name: str, root: Path = PROJECT_ROOT) -> Path:
    path: Path | None = try_find_file_in_project(
        file_name=file_name,
        root=root,
    )

    if path is None:
        raise FileNotFoundError(f"File {file_name!r} was not found under {root.resolve()!s}")

    return path

def resolve_file_from_root(path: str) -> Path:
    resolved_path: Path = PROJECT_ROOT / path
    if not resolved_path.exists():
        raise FileNotFoundError(
            f"File does not exist: {resolved_path}"
        )

    return resolved_path.resolve()

def resolve_dir_from_root(path: str) -> Path:
    resolved_path: Path = PROJECT_ROOT / path

    if not resolved_path.exists():
        raise FileNotFoundError(
            f"Directory does not exist: {resolved_path}"
        )

    if not resolved_path.is_dir():
        raise NotADirectoryError(
            f"Path exists but is not a directory: {resolved_path}"
        )

    return resolved_path.resolve()



def collect_files_in_folder(folder: Path) -> list[Path]:
    if not folder.exists():
        raise FileNotFoundError(f"Folder does not exist: {folder}")

    if not folder.is_dir():
        raise NotADirectoryError(f"Not a directory: {folder}")

    return sorted(
        path
        for path in folder.rglob("*")
        if path.is_file()
    )