import networkx as nx
from typing import Set
from pydantic import Field
from dataclasses import dataclass


@dataclass
class PropagationContainer:
    source_of: Set[str] = Field(default_factory=set)
    # target_of: Set[str] = Field(default_factory=set)  #TODO delete if we officially give up on target sets

class PropagationNetwork(nx.Graph):
    CONTAINER_KEY = "container"
    EDGE_WEIGHT = "weight"

    @classmethod
    def node_has_propagation_container(cls, node_properties):
        return isinstance(node_properties.get(cls.CONTAINER_KEY, None), PropagationContainer)

    def is_ready_for_propagation(self):
        nodes_condition = all(self.node_has_propagation_container(data) for node, data in self.nodes(data=True))
        if not nodes_condition:
            return False
        edges_condition = all(data.get("weight", -1) > 0 for _, _, data in self.edges(data=True))

        return nodes_condition and edges_condition