import math


def log2_transform(values: list[float]) -> list[float]:
    return [
        math.log2(value + 1.0)
        for value in values
    ]