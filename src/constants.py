from enum import Enum

class SpeciesIDs(Enum):
    HUMAN = "human"
    CORONAVIRUS = "sars-cov-2"

class NodeAttrs(Enum):
    SPECIES_ID = "species_id"