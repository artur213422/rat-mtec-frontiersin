import re
from pathlib import Path



def gsm_sample_name_from_cel(path: Path) -> str:
    match = re.match(r"^(GSM\d+)", path.name)

    if match is None:
        raise ValueError(f"No GSM accession found in CEL filename: {path.name}")

    return match.group(1)