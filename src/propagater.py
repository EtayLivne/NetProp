from math import sqrt
from copy import deepcopy
from itertools import chain
from dataclasses import dataclass, field
from functools import total_ordering
import networkx as nx

from datastructures.graph_datastructures import Protein


@dataclass(frozen=False)
class PropagationContainer:
    source_of: set = field(default_factory=set)
    target_of: set = field(default_factory=set)


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


class PropagationError(BaseException):
    pass


class Propagater:
    _PREV_LIQUIDS = "prev_liquids"
    _LIQUIDS = "liquids"
    _EDGE_WEIGHT = "weight"
    _UNNORMALIZED_EDGE_WEIGHT = "unnormalized_weight"
    _DEFAULT_HALT_CONDITION_GAP = 10e-5

    @staticmethod
    def _validate_parameters(propagation_graph, source_condifence_coef, halt_condition_gap):
        if not isinstance(propagation_graph, PropagationGraph):
            raise TypeError(f"propagation_graph must be of type {type(PropagationGraph)}")
        if not (0 < source_condifence_coef < 1):
            raise ValueError("source_confidence_coef must be a number greater than 0 and smaller than 1")
        if not (halt_condition_gap > 0):
            raise ValueError("halt_condition_gap must be a positive number")

    @staticmethod
    def _edge_weight_coef(self, node1_deg, node2_deg):
        return float(1) / (node1_deg + node2_deg)

    def __init__(self, propagation_graph, source_condifence_coef, halt_condition_gap=_DEFAULT_HALT_CONDITION_GAP):
        self._validate_parameters(propagation_graph, source_condifence_coef, halt_condition_gap)

        self.graph = propagation_graph
        self.source_confidence_coef = source_condifence_coef
        self.halt_condition_gap = halt_condition_gap
        self._reset_liquids()
        self._edges_are_normalized = False
        self.unsuppressed_graph = None

    def get_propagation_result_record(self):
        return {node_id: node_data for node_id, node_data in self.graph.nodes(data=True)}

    def validate_propagation_node_sets(self, *nodes):
        missing_nodes = [node for node in chain(*nodes) if node not in self.graph.nodes]
        if missing_nodes:
            print(f'The following nodes are specified in propagation node_set parameters but are not present in '
                  f'the propagation network: {missing_nodes}')
            return False
        return True

    # assumes graph is directed, otherwise each edge will be normalized twice, which is bad!
    def _normalize_edge_weights(self, reverse=False, apply_confidence_coef=True):
        if reverse != self._edges_are_normalized:
            raise PropagationError("May not normalize edges that are already normalized or reverse an unnormalized edge")

        if reverse:
            for node1, node2, data in self.graph.edges(data=True):
                data[self._EDGE_WEIGHT] = data[self._UNNORMALIZED_EDGE_WEIGHT]
            self._edges_are_normalized = False

        else:
            # degree dict for runtime efficiency (access node degree quickly instead of calculating it again for each edge)
            node_degree = {n: self.graph.degree(n) for n in self.graph.nodes()}
            for node_1, node_2, data in self.graph.edges(data=True):
                data[self._UNNORMALIZED_EDGE_WEIGHT] = data[self._EDGE_WEIGHT]

                edge_confidence_coef = float(1 - self.source_confidence_coef) if apply_confidence_coef else 1.0
                data[self._EDGE_WEIGHT] *= edge_confidence_coef * self._edge_weight_coef(node_degree(1), node_degree(2))

            self._edges_are_normalized = True

    def _suppress_nodes(self, suppressed_nodes=None, reverse=False):

        if reverse != bool(self.unsuppressed_graph):
            raise PropagationError("May not attempt to suppress a multiple concurrent node sets, "
                                   "and may not un-suppress a graph that hasn't been suppressed")

        if reverse:
            self.graph = self.unsuppressed_graph
            self.unsuppressed_graph = None

        else:
            suppressed_nodes = suppressed_nodes or set()
            self.unsuppressed_graph = self.graph
            try:
                self.graph = nx.subgraph(self.graph, [node for node in self.graph.nodes if node not in suppressed_nodes])
            except TypeError:
                self.unsuppressed_graph = None
                raise TypeError(f"failed to suppress {suppressed_nodes} because it is not an iterable of node ids")

    def _reset_liquids(self, nodes=None, liquid_types=None):
        # It is possible to reset only a subset of liquids and/or only a subset of nodes. Default is to reset all.
        nodes_to_reset = nodes or self.graph.nodes
        reset_subgraph = nx.subgraph(self.graph, nodes_to_reset)
        for node, data in reset_subgraph.nodes(data=True):
            for liquid_dict in [self._LIQUIDS, self._PREV_LIQUIDS]:
                if liquid_dict not in data:
                    data[liquid_dict] = dict()
            liquids_to_reset = liquid_types if liquid_types else\
                               set(data[self._LIQUIDS].keys()) | set(data[self._PREV_LIQUIDS].keys())
            reset_dict = {liquid_type: 0 for liquid_type in liquids_to_reset}

            data[self._PREV_LIQUIDS].update(reset_dict)
            data[self._LIQUIDS].update(reset_dict)

    def _snapshot_to_prev_liquid(self, nodes, liquid_types=None):
        snapshot_nodes = nodes or self.graph.nodes
        for node in snapshot_nodes:
            data = self.graph.nodes[node]
            snapshot_liquids = liquid_types or data[self._LIQUIDS]
            data[self._PREV_LIQUIDS].update({liquid_type: data[self._LIQUIDS].get(liquid_types, 0)
                                             for liquid_type in snapshot_liquids})

    def _propagate_once(self, subgraph, liquid_types):
        self._snapshot_to_prev_liquid(subgraph.nodes, liquid_types=liquid_types)

        for node_data in subgraph.nodes(data=True):
            node, data = node_data
            for l_t in liquid_types:
                data[self._LIQUIDS][l_t] = \
                    sum(subgraph.nodes[neighbor][self._PREV_LIQUIDS][l_t] * subgraph[node][neighbor][self._EDGE_WEIGHT]
                        for neighbor in subgraph.neighbors(node))
                if l_t in data[self.graph.CONTAINER_PROPERTIES_KEY].source_of:
                    data[self._LIQUIDS][liquid_type] += self.source_confidence_coef

    def _gap_size(self, subgraph, liquid_type):
        return sqrt(sum(pow(data[self._LIQUIDS][liquid_type] - data[self._PREV_LIQUIDS][liquid_type], 2)
                                        for node, data in subgraph.nodes(data=True)))

    def propagate(self, prior_set, suppressed_nodes=None, normalize_flow=False):
        # validate input
        if not self.validate_propagation_node_sets(self.graph, prior_set):
            raise ValueError("cannot propagate because prior_set refers to nodes that do not exist")

        if not self.graph.is_ready_for_propagation():
            raise Exception("cannot propagate because "
                            "graph is not a valid propagation graph")

        self._suppress_nodes(suppressed_nodes=suppressed_nodes)
        self._normalize_edge_weights(apply_confidence_coef=True)

        # Incorporate propagation inputs into propagation graph
        prior_set_subgraph = nx.subgraph(self.graph, prior_set)
        liquid_types = set.union(*[data[self.graph.CONTAINER_PROPERTIES_KEY].source_of for
                                   node, data in prior_set_subgraph.nodes(data=True)])
        self._reset_liquids(liquid_types=liquid_types)

        # Prepare initial state for propagation
        discovered_subgraph_nodes = set(prior_set_subgraph.nodes())
        newest_subgraph_nodes = set(prior_set_subgraph.nodes())
        check_for_gaps_counter = 1
        while True:
            discovered_subgraph = nx.subgraph(self.graph, discovered_subgraph_nodes)
            self._propagate_once(discovered_subgraph, prior_set_subgraph, liquid_types)

            if check_for_gaps_counter % 3 == 0:
                if all(self._gap_size(discovered_subgraph, l_type) < self.halt_condition_gap for l_type in liquid_types):
                    break
            check_for_gaps_counter += 1

            if len(self.graph) > len(discovered_subgraph):
                newest_subgraph_nodes = set.union(*[set(self.graph.neighbors(node)) for node in newest_subgraph_nodes])
                discovered_subgraph_nodes.update(newest_subgraph_nodes)

        # TODO - see if this is actually important/needed, if not: delete.
        # if normalize_flow:
        #     for liquid_type in liquid_types:
        #         type_max = max([node.liquids.get(liquid_type, 0) for node in self.graph.nodes if node not in prior_set])
        #         for node in self.graph.nodes:
        #             node.liquids[liquid_type] = float(node.liquids.get(liquid_type, 0)) / type_max

        # restore suppressed nodes
        self._suppress_nodes(reverse=True)
        self._normalize_edge_weights(reverse=True)

