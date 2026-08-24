from src.common.types.enums import Colors

DATASET_COLORS: dict[str, Colors] = {
    "mouse_gse85_mtec_high": Colors.BLUE,
    "mouse_gse2585_mtec_high": Colors.GREEN,
    "mouse_gse14365_mtec_high": Colors.ORANGE,
    "mouse_gse33878_mtec_high": Colors.PURPLE,
    "rat_mtec_high": Colors.RED,

    "mouse_gse2585_mtec_low":Colors.GREEN,
    "rat_mtec_low": Colors.RED,
    "mouse_gse65617_mtec_low": Colors.ORANGE,
    "mouse_gse205410_mtec_low": Colors.BLUE,
}


def get_dataset_color(dataset_id: str) -> Colors:
    try:
        return DATASET_COLORS[dataset_id]
    except KeyError as error:
        raise ValueError(f"No plot color configured for dataset: {dataset_id}") from error