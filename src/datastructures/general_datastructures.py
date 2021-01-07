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


