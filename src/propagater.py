from math import sqrt
from copy import deepcopy
from itertools import chain
from dataclasses import dataclass, field
from functools import total_ordering
import networkx as nx

from datastructures.graph_datastructures import Protein


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

    def __init__(self, propagation_graph, source_condifence_coef, halt_condition_gap=_DEFAULT_HALT_CONDITION_GAP):
        self._validate_parameters(propagation_graph, source_condifence_coef, halt_condition_gap)

        self.graph = propagation_graph
        self.source_confidence_coef = source_condifence_coef
        self.halt_condition_gap = halt_condition_gap * self.source_confidence_coef
        self._reset_liquids()
        self._edges_are_normalized = False

    @staticmethod
    def _validate_parameters(propagation_graph, source_condifence_coef, halt_condition_gap):
        if not isinstance(propagation_graph, PropagationGraph):
            raise TypeError(f"propagation_graph must be of type {type(PropagationGraph)}")
        if not (0 < source_condifence_coef < 1):
            raise ValueError("source_confidence_coef must be a number greater than 0 and smaller than 1")
        if not (halt_condition_gap > 0):
            raise ValueError("halt_condition_gap must be a positive number")

    # TODO remove once certain that get_propagation_results makes more sense
    # def gene_knockout_copy(self, gene_knockout_ids):
    #     knockout_graph = deepcopy(self.graph)
    #     knockout_graph.remove_nodes_from({Protein(id=gene_id) for gene_id in gene_knockout_ids})
    #     return Propagater(knockout_graph, self.source_confidence_coef, halt_condition_gap=self.halt_condition_gap)

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
        if reverse:
            for node1, node2, data in self.graph.edges(data=True):
                data[self._EDGE_WEIGHT] = data.get(self._UNNORMALIZED_EDGE_WEIGHT, self._EDGE_WEIGHT)
            self._edges_are_normalized = False

        else:
            if self._edges_are_normalized:
                raise PropagationError("May not normalize edges in graph that are already normalized")
            debug_counter = 0
            for node_1, node_2, data in self.graph.edges(data=True):
                norm_factor = sqrt(self.graph.degree(node_1) + self.graph.degree(node_2))
                confidence_coef = (1 - self.source_confidence_coef) if apply_confidence_coef else 1
                data[self._UNNORMALIZED_EDGE_WEIGHT] = data[self._EDGE_WEIGHT]
                data[self._EDGE_WEIGHT] *= confidence_coef / (norm_factor**2)
                if debug_counter % 1000 == 0:
                    print(f'debug counter: {debug_counter}')
                debug_counter += 1
            self._edges_are_normalized = True

    def _get_prev_liquid(self, node, liquid_type):
        return self.graph.nodes[node][self._PREV_LIQUIDS].get(liquid_type, 0)

    def _set_prev_liquid(self, node, prev_liquid, liquid_type):
        self.graph.nodes[node][self._PREV_LIQUIDS][liquid_type] = prev_liquid

    def _suppress_nodes(self, suppressed_nodes=None):
        g = self.graph
        try:
            suppressed_nodes = suppressed_nodes or set()
            self.graph = nx.subgraph(self.graph, [node for node in self.graph.nodes if node not in suppressed_nodes])
        except TypeError:
            raise TypeError("cannot propagate because non-iterable suppressed_nodes provided")
        self._normalize_edge_weights()
        return g

    def _unsuppress_nodes(self, original_graph_pointer):
        self.graph = original_graph_pointer
        self._normalize_edge_weights(reverse=True)

    def _reset_liquids(self, nodes=None, liquid_types=None):
        nodes_to_reset = nodes or self.graph.nodes
        reset_subgraph = nx.subgraph(self.graph, nodes_to_reset)
        for node, data in reset_subgraph.nodes(data=True):
            for liquid_dict in [self._LIQUIDS, self._PREV_LIQUIDS]:
                if liquid_dict not in data:
                    data[liquid_dict] = dict()
            liquids_to_reset = liquid_types if liquid_types else\
                               data[self._LIQUIDS].keys() | data[self._PREV_LIQUIDS].keys()
            reset_dict = {liquid_type: 0 for liquid_type in liquids_to_reset}

            data[self._PREV_LIQUIDS].update(reset_dict)
            data[self._LIQUIDS].update(reset_dict)

    # Stores a snapshot of current liquid amount in each node into the node's prev_liquid property
    def _snapshot_to_prev_liquid(self, node, liquid_types):
        data = self.graph.nodes[node]
        data[self._PREV_LIQUIDS].update({liquid_type: liquid for liquid_type, liquid in data[self._LIQUIDS].items()
                                         if liquid_type in liquid_types})

    def _propagate_once(self, subgraph, prior_set_subgraph, liquid_types):
        for node in subgraph:
            self._snapshot_to_prev_liquid(node, liquid_types=liquid_types)

        for node_data in subgraph.nodes(data=True):
            node, data = node_data
            if not subgraph.degree(node):
                continue
            neighbors = subgraph.neighbors(node)
            data[self._LIQUIDS].\
                update({liquid_type: sum(subgraph.nodes[neighbor][self._PREV_LIQUIDS][liquid_type] * subgraph[node][neighbor][self._EDGE_WEIGHT]
                                         for neighbor in neighbors)
                        for liquid_type in liquid_types})

        # Refill prior set
        for node, data in prior_set_subgraph.nodes(data=True):
            for liquid_type in data[self.graph.CONTAINER_PROPERTIES_KEY].source_of:
                data[self._LIQUIDS][liquid_type] += self.source_confidence_coef

    def propagate(self, prior_set, suppressed_nodes=None, normalize_flow=False):
        # validate input
        if not self.validate_propagation_node_sets(self.graph, prior_set):
            raise ValueError("cannot propagate because prior_set refers to nodes that do not exist")

        if not self.graph.is_ready_for_propagation():
            raise Exception("cannot propagate because "
                            "graph is not a valid propagation graph")

        g = self._suppress_nodes(suppressed_nodes=suppressed_nodes)

        # Incorporate propagation inputs into propagation graph
        prior_set_subgraph = nx.subgraph(self.graph, prior_set)
        liquid_types = set.union(*[data[self.graph.CONTAINER_PROPERTIES_KEY].source_of for
                                   node, data in prior_set_subgraph.nodes(data=True)])
        self._reset_liquids(liquid_types=liquid_types)

        # Prepare initial state for propagation
        discovered_subgraph_nodes = set(prior_set_subgraph.nodes())
        newest_subgraph_nodes = set(prior_set_subgraph.nodes())
        gaps = {liquid_type: 10*self.halt_condition_gap for liquid_type in liquid_types}
        check_for_gaps_counter = 1
        while True:
            print("number of prop iterations: {}".format(check_for_gaps_counter))
            discovered_subgraph = nx.subgraph(self.graph, discovered_subgraph_nodes)
            self._propagate_once(discovered_subgraph, prior_set_subgraph, liquid_types)

            if check_for_gaps_counter % 3 == 0:
                for l_type in liquid_types:
                    gaps[l_type] = sqrt(sum(pow(data[self._LIQUIDS][l_type] - data[self._PREV_LIQUIDS][l_type], 2)
                                        for node, data in discovered_subgraph.nodes(data=True)))
                if max(gaps.values()) < self.halt_condition_gap:
                    break
                print(f'smallest gap is: {min(gaps.values())}')
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
        self._unsuppress_nodes(g)

