from enum import StrEnum


class Genotypes(StrEnum):
    CTRL = "ctrl"
    KO  = "ko"



class Species(StrEnum):
    RAT = "rat"
    MOUSE = "mouse"

class CellPopulation(StrEnum):
    MTEC_HIGH = "mTEC_high"
    MTEC_LOW = "mTEC_low"
