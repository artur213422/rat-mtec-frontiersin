from enum import StrEnum


class Genotypes(StrEnum):
    CTRL = "ctrl"
    KO  = "ko"
    ALL = "all"



class Species(StrEnum):
    RAT = "rat"
    MOUSE = "mouse"

class CellPopulation(StrEnum):
    MTEC_HIGH = "mTEC_high"
    MTEC_LOW = "mTEC_low"


class Colors(StrEnum):
    RED = "red"
    GREEN = "green"
    BLUE = "blue"
    BLACK = "black"
    ORANGE = "orange"
    PURPLE = "purple"
    WHITE = "white"