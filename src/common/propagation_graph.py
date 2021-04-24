from dataclasses import dataclass, field
from functools import total_ordering
import networkx as nx


@dataclass(frozen=False)
@total_ordering
class PropagationContainer:
    id: int
    source_of: set = field(default_factory=set)
    target_of: set = field(default_factory=set)

    def __hash__(self):
        return self.id

    def __eq__(self, other):
        return self.id == other.id

    def __lt__(self, other):
        return self.id < other.id


class PropagationResult(PropagationContainer):
    liquids: set = field(default_factory=set)
    in_prior_set: bool


class PropagationGraph(nx.Graph):
    CONTAINER_PROPERTIES_KEY = "container_properties"

    @classmethod
    def node_has_propagation_container(cls, node_properties):
        return isinstance(node_properties.get(cls.CONTAINER_PROPERTIES_KEY, None), PropagationContainer)

    def is_ready_for_propagation(self):
        nodes_condition = all(self.node_has_propagation_container(data) for node, data in self.nodes(data=True))
        if not nodes_condition:
            return False
        ids_condition = all(node == data[self.CONTAINER_PROPERTIES_KEY].id for node, data in self.nodes(data=True))
        edges_condition = all(data.get("weight", -1) > 0 for source, dest, data in self.edges(data=True))

        return nodes_condition and ids_condition and edges_condition
