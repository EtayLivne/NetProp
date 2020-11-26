from dataclasses import dataclass, field

@dataclass()
class CovTargetMetadata:
    name: str
    sources: set = field(default_factory=set)
    infection_roles: set = field(default_factory=set)

    def add_sources(self, sources):
        if type(sources) is str:
            self.sources.update({sources})
        else:
            self.sources.update(sources)

    def add_infection_roles(self, roles):
        if type(roles) is str:
            self.infection_roles.update({roles})
        else:
            self.infection_roles.update(roles)

