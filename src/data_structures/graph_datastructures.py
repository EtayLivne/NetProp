from dataclasses import dataclass
from functools import total_ordering


@dataclass(frozen=False)
@total_ordering
class Protein:
    id: int
    liquid: float = 0
    rank: int = 0

    def __hash__(self):
        return self.id

    def __eq__(self, other):
        return self.id == other.id

    def __lt__(self, other):
        return self.liquid < other.liquid
