from math import sqrt
from itertools import chain
from pydantic.dataclasses import dataclass
from pydantic import Field
from typing import Set
#from dataclasses import dataclass, field
import networkx as nx


@dataclass
class PropagationContainer:
    source_of: Set[str] = Field(default_factory=set)
    target_of: Set[str] = Field(default_factory=set)


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


class PropagationError(BaseException):
    pass

class PropagationNetworkError(BaseException):
    pass

class Propagater:
    _PREV_LIQUIDS = "prev_liquids"
    _LIQUIDS = "liquids"
    _EDGE_WEIGHT = PropagationNetwork.EDGE_WEIGHT
    _UNNORMALIZED_EDGE_WEIGHT = "unnormalized_weight"
    _DEFAULT_HALT_CONDITION_GAP = 10e-5
    NO_MAX_STEPS = -1
    NO_MIN_GAP = -1

    @staticmethod
    def _edge_degree_coef(node1_deg, node2_deg):
        return float(2) / (node1_deg + node2_deg)

    def __init__(self, network, source_confidence, min_gap=NO_MIN_GAP, max_steps=NO_MAX_STEPS):
        self.network = network
        self.source_confidence = source_confidence
        self.max_steps = max_steps
        self.min_gap = min_gap

        self._reset_liquids()
        self._edges_are_normalized = False
        self._unsuppressed_network = None

    @property
    def network(self):
        return self._network

    @network.setter
    def network(self, network):
        if not isinstance(network, PropagationNetwork):
            raise TypeError(f"Propagation_network must be of type {type(PropagationNetwork)}, given {type(network)}")
        if not network.is_ready_for_propagation():
            raise PropagationNetworkError("Network does not conform to all structural requirements for propagation")
        self._network = network

    @property
    def source_confidence(self):
        return self._source_confidence

    @source_confidence.setter
    def source_confidence(self, source_confidence):
        if not (0 < source_confidence < 1):
            raise ValueError(f"source_confidence must be a real number between 1 and 0, given {source_confidence}")
        self._source_confidence = source_confidence

    @property
    def max_steps(self):
        return self._max_steps

    @max_steps.setter
    def max_steps(self, max_steps):
        if max_steps <= 0 and max_steps != self.NO_MAX_STEPS:
            raise ValueError(f"if max number of steps is defined, it must be a positive integer, given {max_steps}")
        self._max_steps = max_steps

    @property
    def min_gap(self):
        return self._min_gap

    @min_gap.setter
    def min_gap(self, min_gap):
        if min_gap <= 0 and min_gap != self.NO_MIN_GAP:
            raise ValueError(f"min_gap must be a positive number, given {min_gap}")
        self._min_gap = min_gap

    def validate_propagation_node_sets(self, *nodes):
        missing_nodes = [node for node in chain(*nodes) if node not in self._network.nodes]
        if missing_nodes:
            print(f'The following nodes are specified in propagation node_set parameters but are not present in '
                  f'the propagation network: {missing_nodes}')
            return False
        return True

    def propagate(self, prior_set, suppressed_nodes=None):
        # validate input
        if not self.validate_propagation_node_sets(self._network, prior_set):
            raise ValueError("Cannot propagate because prior_set refers to nodes that do not exist")

        if self._max_steps == self.NO_MAX_STEPS and self._min_gap == self.NO_MIN_GAP:
            raise PropagationError("No halt condition specified")

        # Incorporate propagation inputs into propagation graph
        self._suppress_nodes(suppressed_nodes=suppressed_nodes)
        self._normalize_edge_weights(apply_confidence_coef=True)
        liquid_types = set.union(*[self._network.nodes[node][self._network.CONTAINER_KEY].source_of
                                   for node in prior_set if node in prior_set])
        self._reset_liquids(liquid_types=liquid_types)

        # Prepare initial propagation state
        discovered_subgraph_nodes = prior_set
        newest_subgraph_nodes = discovered_subgraph_nodes
        check_for_gaps_counter = 1
        step_counter = 1

        # Propagate
        while True:
            discovered_subgraph = nx.subgraph(self._network, discovered_subgraph_nodes)
            if len(self._network) > len(discovered_subgraph):
                newest_subgraph_nodes = set.union(*[set(self._network.neighbors(node)) for node in newest_subgraph_nodes])
                discovered_subgraph_nodes.update(newest_subgraph_nodes)

            self._propagate_once(discovered_subgraph, prior_set, liquid_types)

            # Test halt conditions
            if self._max_steps != self.NO_MAX_STEPS:
                if step_counter >= self.max_steps:
                    break
                step_counter += 1
            elif self._min_gap != self.NO_MIN_GAP and check_for_gaps_counter % 3 == 0:
                if all(self._gap_size(discovered_subgraph, l_type) < self._min_gap for l_type in liquid_types):
                    break
                check_for_gaps_counter += 1

        # restore suppressed nodes
        self._suppress_nodes(reverse=True)
        self._normalize_edge_weights(reverse=True)

    def node_liquids(self, node_id):
        return self.network[node_id][self._LIQUIDS]

    def _reset_liquids(self, nodes=None, liquid_types=None):
        # It is possible to reset only a subset of liquids and/or only a subset of nodes. Default is to reset all.
        nodes_to_reset = nodes or self._network.nodes
        reset_subgraph = nx.subgraph(self._network, nodes_to_reset)
        for node, data in reset_subgraph.nodes(data=True):
            for liquid_dict in [self._LIQUIDS, self._PREV_LIQUIDS]:
                if liquid_dict not in data:
                    data[liquid_dict] = dict()
            liquids_to_reset = liquid_types if liquid_types else\
                               set(data[self._LIQUIDS].keys()) | set(data[self._PREV_LIQUIDS].keys())
            reset_dict = {liquid_type: 0 for liquid_type in liquids_to_reset}

            data[self._PREV_LIQUIDS].update(reset_dict)
            data[self._LIQUIDS].update(reset_dict)

    # assumes graph is directed, otherwise each edge will be normalized twice, which is bad!
    def _normalize_edge_weights(self, reverse=False, apply_confidence_coef=True):
        if reverse != self._edges_are_normalized:
            raise PropagationError(
                "May not normalize edges that are already normalized or reverse an unnormalized edge")

        if reverse:
            for node1, node2, data in self._network.edges(data=True):
                data[self._EDGE_WEIGHT] = data[self._UNNORMALIZED_EDGE_WEIGHT]
            self._edges_are_normalized = False

        else:
            # degree dict to access node degree quickly instead of calculating it again for each edge
            node_degree = {n: self._network.degree(n) for n in self._network.nodes()}
            confidence_coef = float(1 - self._source_confidence) if apply_confidence_coef else 1.0
            for node_1, node_2, data in self._network.edges(data=True):
                data[self._UNNORMALIZED_EDGE_WEIGHT] = data[self._EDGE_WEIGHT]
                data[self._EDGE_WEIGHT] *=\
                    confidence_coef * self._edge_degree_coef(node_degree[node_1], node_degree[node_2])

            self._edges_are_normalized = True

    def _suppress_nodes(self, suppressed_nodes=None, reverse=False):

        if reverse and not bool(self._unsuppressed_network):
            raise PropagationError("No nodes to un-suppress")
        if bool(self._unsuppressed_network) and not reverse:
            raise PropagationError("May not suppress nodes while another set of nodes is already supressed")

        if reverse:
            self._network = self._unsuppressed_network
            self._unsuppressed_network = None

        else:
            suppressed_nodes = suppressed_nodes or set()
            suppressed_nodes_list = [node for node in self._network.nodes if node not in suppressed_nodes]

            self._unsuppressed_network = self._network
            try:
                self._network = nx.subgraph(self._network, suppressed_nodes_list)
            except TypeError:
                self._unsuppressed_network = None
                raise TypeError(f"failed to suppress {suppressed_nodes} because it is not an iterable of node ids")

    def _propagate_once(self, subgraph, prior_set, liquid_types):
        self._snapshot_to_prev_liquid(subgraph.nodes, liquid_types=liquid_types)

        for node_data in subgraph.nodes(data=True):
            node, data = node_data
            for l_t in liquid_types:
                data[self._LIQUIDS][l_t] = \
                    sum(subgraph.nodes[neighbor][self._PREV_LIQUIDS][l_t] * subgraph[node][neighbor][self._EDGE_WEIGHT]
                        for neighbor in subgraph.neighbors(node))
                if node in prior_set and l_t in data[self._network.CONTAINER_KEY].source_of:
                    data[self._LIQUIDS][l_t] += self._source_confidence

    def _snapshot_to_prev_liquid(self, nodes, liquid_types=None):
        snapshot_nodes = nodes or self._network.nodes
        for node in snapshot_nodes:
            data = self._network.nodes[node]
            snapshot_liquids = liquid_types or data[self._LIQUIDS]
            data[self._PREV_LIQUIDS].update({liquid_type: data[self._LIQUIDS].get(liquid_type, 0)
                                             for liquid_type in snapshot_liquids})

    def _gap_size(self, subgraph, liquid_type):
        return sqrt(sum(pow(data[self._LIQUIDS][liquid_type] - data[self._PREV_LIQUIDS][liquid_type], 2)
                                        for node, data in subgraph.nodes(data=True)))


network = PropagationNetwork()
network.add_edge("1", "2")
network.nodes["1"][network.CONTAINER_KEY] = PropagationContainer(source_of={"liquid"})
network.nodes["2"][network.CONTAINER_KEY] = PropagationContainer(target_of={"liquid"})
network["1"]["2"][network.EDGE_WEIGHT] = 1
p = Propagater(network, 0.5, max_steps=3)
p.propagate({"1"})