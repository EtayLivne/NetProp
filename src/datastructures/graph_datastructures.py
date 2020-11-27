from dataclasses import dataclass, field
from functools import total_ordering


@dataclass(frozen=False)
@total_ordering
class Protein:
    id: int
    liquids: dict = field(default_factory=dict)
    source_of: set = field(default_factory=set)
    rank: int = 0

    def __hash__(self):
        return self.id

    def __eq__(self, other):
        return self.id == other.id

    def __lt__(self, other):
        return self.liquid < other.liquid


def construct_prior_set(prior_set_seed):
    if type(prior_set_seed) is dict:
        try:
            return {Protein(id=int(k), source_of=set(v)) for k, v in prior_set_seed.items()}
        except ValueError:
            raise ValueError("keys of prior set seed dict must be integers and values must be iterables")
    else:
        try:
            return {Protein(id=int(prior), source_of={int(prior)}) for prior in iter(prior_set_seed)}
        except ValueError:
            raise ValueError("members of iterable prior set seed must be integers")
        except TypeError:
            raise TypeError("prior set seed must be a dict or an iterable!")