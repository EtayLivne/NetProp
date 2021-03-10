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


class PropagationGraph(nx.DiGraph):
    CONTAINER_PROPERTIES_KEY = "container_properties"

    @classmethod
    def node_has_propagation_container(cls, node_properties):
        return isinstance(node_properties.get(cls.CONTAINER_PROPERTIES_KEY, None), PropagationContainer)

    def is_ready_for_propagation(self):
        node_tuple_id_index = 0
        node_tuple_properties_index = 1
        edge_tuple_properties_index = 2
        nodes_condition = all(self.node_has_propagation_container(node_tuple[node_tuple_properties_index])
                              for node_tuple in self.nodes(data=True))
        ids_condition = all(node_tuple[node_tuple_id_index] == node_tuple[node_tuple_properties_index].get("id", None)
                            for node_tuple in self.nodes(data=True))
        edges_condition = all(edge_tuple[edge_tuple_properties_index].get("weight", -1) > 0
                              for edge_tuple in self.edges(data=True))

        return nodes_condition and ids_condition and edges_condition



class Propagater:
    _PREV_LIQUIDS_PROPERTY_KEY = "prev_liquids"
    _DEFAULT_HALT_CONDITION_GAP = 10e-2
    _EDGE_WEIGHT_KEY = "weight"
    _UNNORMALIZED_EDGE_WEIGHT_KEY = "unnormalized_weight"

    def __init__(self, propagation_graph, source_condifence_coef, halt_condition_gap=_DEFAULT_HALT_CONDITION_GAP):
        self._validate_parameters(propagation_graph, source_condifence_coef, halt_condition_gap)

        self.graph = propagation_graph
        self.source_confidence_coef = source_condifence_coef
        self.halt_condition_gap = halt_condition_gap
        self._reset_liquids()

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
        missing_nodes = [node.id for node in chain(*nodes) if node not in self.graph.nodes.keys()]
        if missing_nodes:
            print(f'The following nodes are specified in propagation node_set parameters but are not present in '
                  f'the propagation network: {missing_nodes}')
            return False
        return True

    # assumes graph is directed, otherwise each edge will be normalized twice, which is bad!
    def _normalize_edge_weights(self, reverse=False):
        if reverse:
            for edge_data in self.graph.edges(data=True):
                edge_data[2][self._EDGE_WEIGHT_KEY] = edge_data[2][self._UNNORMALIZED_EDGE_WEIGHT_KEY]
        else:
            for node in self.graph:
                edges = [self.graph[node][neighbor] for neighbor in self.graph.neighbors(node)]
                sum_of_weights = sum(edge[self._EDGE_WEIGHT_KEY] for edge in edges)
                for edge in edges:
                    edge[self._UNNORMALIZED_EDGE_WEIGHT_KEY] = edge[self._EDGE_WEIGHT_KEY]
                    edge[self._EDGE_WEIGHT_KEY] = float(edge[self._EDGE_WEIGHT_KEY]) / sum_of_weights

    def _get_prev_liquid(self, node, liquid_type):
        return self.graph.nodes[node][self._PREV_LIQUIDS_PROPERTY_KEY].get(liquid_type, 0)

    def _set_prev_liquid(self, node, prev_liquid, liquid_type):
        self.graph.nodes[node][self._PREV_LIQUIDS_PROPERTY_KEY][liquid_type] = prev_liquid

    def _reset_liquids(self, nodes=None, liquid_types=None):
        nodes = nodes if nodes else self.graph.nodes
        for node in nodes:
            liquids_to_reset = liquid_types if liquid_types else node.liquids.keys()
            reset_dict = {liquid_type: 0 for liquid_type in liquids_to_reset}

            self.graph.nodes[node][self._PREV_LIQUIDS_PROPERTY_KEY].update(reset_dict)
            node.liquids.update(reset_dict)

    # Stores a snapshot of current liquid amount in each node into the node's prev_liquid property
    def _snapshot_to_prev_liquid(self, node, liquid_types):
        self.graph.nodes[node][self._PREV_LIQUIDS_PROPERTY_KEY].\
            update({liquid_type: liquid for liquid_type, liquid in node.liquids if liquid_type in liquid_types})

    def _propagate_once(self, subgraph, prior_set_subgraph, liquid_types):
        alpha = self.source_confidence_coef
        for node in subgraph:
            self._snapshot_to_prev_liquid(node, liquid_types=liquid_types)

        # If liquids dict is empty in any node in prior set then this is the first propagation step
        for node in subgraph.nodes():
            for neighbor in subgraph.neighbors(node):
                edge_capacity = subgraph[node][neighbor][self._EDGE_WEIGHT_KEY] * alpha
                for liquid_type in liquid_types:
                    neighbor.liquids[liquid_type] = neighbor.liquids.get(liquid_type, 0) + \
                                                    self._get_prev_liquid(node, liquid_type) * edge_capacity
            for liquid_type in liquid_types:
                node.liquids[liquid_type] = node.liquids.get(liquid_type, 0) - self._get_prev_liquid(node, liquid_type)
            # pump new liquid for next round

        for node in prior_set_subgraph:
            node.liquids.update({liquid_type: node.liquids[liquid_type] + 1 - alpha for liquid_type in node.source_of})

    def propagate(self, prior_set, suppressed_nodes=None, normalize_flow=False):
        # validate input
        if not self.validate_propagation_node_sets(self.graph, prior_set):
            raise ValueError("cannot propagate because prior_set refers to nodes that do not exist")

        if not self.graph.is_ready_for_propagation():
            raise Exception("cannot propagate because "
                            "graph is not a valid propagation graph")

        # suppress nodes
        g = self.graph
        try:
            suppressed_nodes = suppressed_nodes or set()
            self.graph = nx.subgraph(self.graph, [node for node in self.graph.nodes if node not in suppressed_nodes])
        except TypeError:
            raise TypeError("cannot propagate because non-iterable suppressed_nodes provided")
        self._normalize_edge_weights()

        # Incorporate propagation inputs into propagation graph
        liquid_types = set.union(*[prior.source_of for prior in prior_set])
        self._reset_liquids(liquid_types=liquid_types)
        prior_set_subgraph = nx.subgraph(self.graph, [node.id for node in self.graph.nodes.keys()
                                                      if node in prior_set and node.id not in suppressed_nodes])

        # Prepare initial state for propagation
        discovered_subgraph_nodes = set(prior_set_subgraph.nodes())
        newest_subgraph_nodes = set(prior_set_subgraph.nodes())
        gaps = {liquid_type: 10*self.halt_condition_gap for liquid_type in liquid_types}
        check_for_gaps_counter = 1
        while True:
            discovered_subgraph = nx.subgraph(self.graph, discovered_subgraph_nodes)
            self._propagate_once(discovered_subgraph, prior_set_subgraph, liquid_types)

            if check_for_gaps_counter % 3 == 0:
                for l_type in liquid_types:
                    gaps[l_type] = sqrt(sum(pow(node.liquids[l_type] - self._get_prev_liquid(node, l_type), 2)
                                        for node in discovered_subgraph.nodes.keys()))
                if min(gaps.values()) < self.halt_condition_gap:
                    break
                print(f'smallest gap is: {min(gaps.values())}')
            else:
                check_for_gaps_counter += 1

            if len(self.graph) > len(discovered_subgraph):
                newest_subgraph_nodes = set.union(*[set(self.graph.neighbors(node)) for node in newest_subgraph_nodes])
                discovered_subgraph_nodes.update(newest_subgraph_nodes)

        if normalize_flow:
            for liquid_type in liquid_types:
                type_max = max([node.liquids.get(liquid_type, 0) for node in self.graph.nodes if node not in prior_set])
                for node in self.graph.nodes:
                    node.liquids[liquid_type] = float(node.liquids.get(liquid_type, 0)) / type_max

        # restore suppressed nodes
        self._normalize_edge_weights(reverse=True)
        self.graph = g
