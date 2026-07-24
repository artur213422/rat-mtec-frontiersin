from src.common.types.enums import Colors

DATASET_COLORS: dict[str, Colors] = {
    "mouse_gse85": Colors.BLUE,
    "mouse_gse2585": Colors.GREEN,
    "mouse_gse14365": Colors.ORANGE,
    "mouse_gse33878": Colors.PURPLE,
    "rat_mtec_high": Colors.RED,
}


def get_dataset_color(dataset_id: str) -> Colors:
    try:
        return DATASET_COLORS[dataset_id]
    except KeyError as error:
        raise ValueError(
            f"No plot color configured for dataset: {dataset_id}"
        ) from error