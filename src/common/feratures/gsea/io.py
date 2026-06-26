from pathlib import Path


def gmt_to_gene_sets_dict(gmt_file: Path) -> dict[str, list[str]]:
    gene_sets: dict[str, list[str]] = {}

    with gmt_file.open("r", encoding="utf-8") as file:
        for line_number, line in enumerate(file, start=1):
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 3:
                raise ValueError(
                    f"Invalid GMT file '{gmt_file}', line {line_number}: "
                    f"expected at least 3 tab-separated columns."
                )

            pathway: str = fields[0]
            genes: list[str] = sorted(set(fields[2:]))

            gene_sets[pathway] = genes

    return gene_sets