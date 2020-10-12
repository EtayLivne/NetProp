from dataclasses import dataclass
import networkx as nx

@dataclass
class Protein:
    id: int
    liquid: float
    rank: int