from enum import Enum


class Namespace(Enum):
    BP = "biological_process"
    MF = "molecular_function"
    CP = "cellular_component"


def to_namespace(namespace: str) -> Namespace:
    try:
        return Namespace(namespace)
    except ValueError as error:
        raise ValueError(f"Unknown GO namespace: '{namespace}'") from error


