
def _get_legend_label(label: str | None, used_labels: set[str]) -> str | None:
    if label is None:
        return None

    if label in used_labels:
        return None

    used_labels.add(label)

    return label