import json
from dataclasses import dataclass, field

@dataclass()
class PriorMetadata:
    """
    represents metadata about a single member of the prior set for ppi network propagation purposes.
    name: the standard symbol representing this protein. If cov protein, the symbol is taken from the paper cited in readme.txt
    sources: list of symbols for sources in cov ppi for this prior set member. If cov protein, will usually just contain itself
    infection_roles: set of roles in infection process that this prior member is a part of.
    """
    name: str
    infection_roles: set = field(default_factory=set)

    def add_infection_roles(self, roles):
        if type(roles) is str:
            self.infection_roles.update({roles})
        else:
            self.infection_roles.update(roles)


class HumanPriorMetadata(PriorMetadata):
    sources: set = field(default_factory=set)

    def add_sources(self, sources):
        if type(sources) is str:
            self.sources.update({sources})
        else:
            self.sources.update(sources)


class CovPriorMetadata(PriorMetadata):
    targets: set = field(default_factory=set)

    def add_targets(self, targets):
        if type(targets) is str:
            self.sources.update({targets})
        else:
            self.sources.update(targets)

@dataclass(frozen=True)
class KnockoutGeneSet:
    name: str
    gene_set: set = field(default_factory=set)

    def __iter__(self):
        return iter(self.gene_set)


class Singleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class SymbolEntrezgeneMap(metaclass=Singleton):
    def __init__(self, from_file=None):
        self._symbol_to_entrezgene_map = dict()
        if from_file:
            self._init_from_file(from_file)

    def _init_from_file(self, path):
        with open(path, 'r') as json_handler:
            self._symbol_to_entrezgene_map = json.load(json_handler)

    def get_entrezgene(self, symbol):
        return self._symbol_to_entrezgene_map[symbol]

    def get_symbol(self, entrezgene):
        try:
            return next(symbol for symbol in self._symbol_to_entrezgene_map if self._symbol_to_entrezgene_map[symbol] == entrezgene)
        except StopIteration:
            raise KeyError


def init_symbol_entrezgene_map(from_file=None):
    return SymbolEntrezgeneMap(from_file=from_file)