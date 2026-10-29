from pathlib import Path


# Folder containing files to export
input_folder: Path = Path(__file__).parent

# Output markdown file
output_file: Path = input_folder / "manuscript.md"

# File types worth sending to AI
INCLUDE_EXTENSIONS = {
    ".py",
    ".cpp",
    ".hpp",
    ".h",
    ".c",
    ".cc",
    ".md",
    ".txt",
    ".json",
    ".yaml",
    ".yml",
    ".toml",
    ".ini",
}


def should_include(path: Path) -> bool:
    # Skip this script itself
    if path.resolve() == Path(__file__).resolve():
        return False

    # Skip generated output
    if path.resolve() == output_file.resolve():
        return False

    # Skip hidden files/folders
    if any(part.startswith(".") for part in path.parts):
        return False

    # Only selected extensions
    return path.suffix.lower() in INCLUDE_EXTENSIONS


def export_context():
    files = sorted(
        p for p in input_folder.rglob("*")
        if p.is_file() and should_include(p)
    )

    with output_file.open("w", encoding="utf-8") as md:
        md.write("# Project Context Export\n\n")
        md.write(f"Source folder: `{input_folder}`\n\n")
        md.write("---\n\n")

        for file in files:
            relative = file.relative_to(input_folder)

            md.write(f"## File: `{relative}`\n\n")

            try:
                content = file.read_text(encoding="utf-8")
            except UnicodeDecodeError:
                md.write("*Binary or unsupported encoding skipped.*\n\n")
                continue

            extension = file.suffix.replace(".", "") or "text"

            md.write(f"```{extension}\n")
            md.write(content)
            if not content.endswith("\n"):
                md.write("\n")
            md.write("```\n\n")
            md.write("---\n\n")

    print(f"Created: {output_file}")
    print(f"Included files: {len(files)}")


if __name__ == "__main__":
    export_context()