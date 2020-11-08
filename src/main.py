import networkx as nx
from copy import deepcopy
from math import sqrt, pow, floor
from collections.abc import Iterable
from itertools import chain

# TODO finish refactor of collection of functions to class
# TODO maybe save both normalized and unnormalized edge weights?
# TODO support multiprocessing, shouldn't be difficult (famous last words, and this is really more of a bonus)
WEIGHT = "weight"
DATA = 2
from data_structures.graph_datastructures import Protein


# TODO optimize away the .get in liquid spilling by initializing the gap keys for all nodes with 0 value for prior set nodes

def graph_from_file(file_path):
    source_node_index = 0
    target_node_index = 1
    edge_weight_index = 2

    graph = nx.DiGraph()
    discovered_nodes_map = dict()
    with open(file_path, 'r') as handler:
        for line in handler.readlines():
            values = line.split()
            for node_index in (values[source_node_index], values[target_node_index]):
                if node_index not in discovered_nodes_map:
                    discovered_nodes_map[node_index] = Protein(id=int(node_index))

            graph.add_edge(discovered_nodes_map[values[source_node_index]],
                           discovered_nodes_map[values[target_node_index]],
                           weight=float(values[edge_weight_index]))
    return graph


def get_propagater(source, from_file=True):
    graph = graph_from_file(source) if from_file else deepcopy(source)
    return Propagater(graph)


class Propagater:
    _PREV_LIQUID_PROPERTY_KEY = "prev_liquids"
    _DEFUALT_HALT_CONDITION_GAP = 10e-3

    def __init__(self, graph, confidence_coefficient, halt_condition_gap=_DEFUALT_HALT_CONDITION_GAP):
        self.graph = graph
        self.confidence_coefficient = confidence_coefficient
        self.halt_condition_gap = halt_condition_gap
        self.gaps = None

        self._normalize_edge_weights()
        self._reset_propagation_parameters()

    @staticmethod
    def validate_propagation_node_sets(network, *node_id_sets):
        missing_nodes = [n_id for n_id in chain(*node_id_sets) if n_id not in [node.id for node in network.nodes()]]
        if missing_nodes:
            print(f'The following nodes are specified in propagation node_set parameters but are not present in '
                  f'the propagation network: {missing_nodes}')

    def _reset_propagation_parameters(self):
        self.gaps = {node: 0 for node in self.graph}
        for node in self.graph:
            node.liquid = 0
            self._reset_liquid(node)

    # assumes graph is directed, otherwise each edge will be normalized twice, which is bad!
    def _normalize_edge_weights(self):
        for node in self.graph:
            edges = [self.graph[node][neighbor] for neighbor in self.graph.neighbors(node)]
            sum_of_weights = sum(edge[WEIGHT] for edge in edges)
            for edge in edges:
                edge[WEIGHT] = float(edge[WEIGHT]) / sum_of_weights

    def __get_node_properties(self, node):
        return self.graph.nodes[node]

    def _get_prev_liquid(self, node, liquid_source):
        node_properties = self.__get_node_properties(node)
        return node_properties[self._PREV_LIQUID_PROPERTY_KEY].get(liquid_source, 0)

    def _set_prev_liquid(self, node, prev_liquid, liquid_source):
        node_properties = self.__get_node_properties(node)
        node_properties[self._PREV_LIQUID_PROPERTY_KEY][liquid_source] = prev_liquid

    def _reset_liquid(self, node, liquid_sources=None):
        node_properties = self.__get_node_properties(node)
        if liquid_sources:
            for liquid_source in liquid_sources:
                del(node_properties[self._PREV_LIQUID_PROPERTY_KEY][liquid_source])
        else:
            node_properties[self._PREV_LIQUID_PROPERTY_KEY] = dict()

    # Stores a snapshot of current liquid amount in each node into the node's prev_liquid property
    def _snapshot_to_prev_liquid(self, node):
        node_prev_liquids = self.__get_node_properties(node)[self._PREV_LIQUID_PROPERTY_KEY]
        node_prev_liquids.update(node.liquids)

    def _propagate_once(self, subgraph, prior_set_subgraph, confidence_coefficient):
        alfa = confidence_coefficient
        for node in subgraph:
            self._snapshot_to_prev_liquid(node)

        # If liquids dict is empty in any node in prior set then this is the first propagation step
        if not next(iter(prior_set_subgraph)).liquids:
            for node in prior_set_subgraph:
                node.liquids[node] = 1 - alfa
                self.gaps[node] = 1 - alfa
        else:
            for node in subgraph.nodes():
                for neighbor in subgraph.neighbors(node):
                    edge_capacity = subgraph[node][neighbor][WEIGHT] * alfa
                    for source_node in prior_set_subgraph:
                        neighbor.liquids[source_node] = neighbor.liquids.get(source_node, 0) +\
                                                        node.liquids.get(source_node, 0) * edge_capacity
                for source_node in prior_set_subgraph:
                    node.liquids[source_node] = node.liquids.get(source_node, 0) - self._get_prev_liquid(node, source_node)
                # pump new liquid for next round
                if node in prior_set_subgraph:
                    node.liquids[node] = node.liquids[node] + (1 - alfa)

            for source_node in prior_set_subgraph:
                self.gaps[source_node] = sqrt(sum(pow(node.liquids[source_node] - self._get_prev_liquid(node, source_node), 2)
                                                  for node in subgraph.nodes()))

    def propagate(self, prior_set_ids):
        self.validate_propagation_node_sets(self.graph, prior_set_ids)
        self._reset_propagation_parameters()
        prior_set = {node for node in self.graph.nodes if node.id in prior_set_ids}
        prior_set_subgraph = nx.subgraph(self.graph, prior_set)

        discovered_subgraph_nodes = set(prior_set)
        newest_subgraph_nodes = set(prior_set)
        # set gaps for all nodes in prior set to >> halt condition
        self.gaps.update({node: 10*self.halt_condition_gap for node in prior_set})
        while [v for v in self.gaps.values() if v > self.halt_condition_gap]:
            discovered_subgraph = nx.subgraph(self.graph, discovered_subgraph_nodes)
            self._propagate_once(discovered_subgraph, prior_set_subgraph, self.confidence_coefficient)
            newest_subgraph_nodes = set.union(*[set(self.graph.neighbors(node)) for node in newest_subgraph_nodes])
            discovered_subgraph_nodes.update(newest_subgraph_nodes)

    # deprecated
    def top_k_candidates(self, k):
        raise NotImplemented
        top_k = floor(len(self.graph) / k)
        return sorted(self.graph.nodes(), reverse=True)[:top_k]


def measure_gene_knockout_impact(network_source_file_path, prior_set_ids, target_set_ids, knockout_candidate_ids):
    if type(knockout_candidate_ids) is int:
        knockout_candidate_ids = {knockout_candidate_ids}
    elif not isinstance(knockout_candidate_ids, Iterable):
        raise TypeError("knockoout_candidate_ids must be a single int or iterable of ints")
    network = graph_from_file(network_source_file_path)
    knockout_network = deepcopy(network)
    knockout_network.remove_nodes_from(knockout_candidate_ids)
    Propagater.validate_propagation_node_sets(network, prior_set_ids, target_set_ids, knockout_candidate_ids)

    confidence_coefficient = 0.3
    network_propagater = Propagater(network, confidence_coefficient)
    knockout_network_propagater = Propagater(knockout_network, confidence_coefficient)

    network_propagater.propagate(prior_set_ids)
    knockout_network_propagater.propagate(prior_set_ids.difference(knockout_candidate_ids))
    # TODO plot differences in graphs



if __name__ == "__main__":
    ppi = graph_from_file(r"..\data\H_sapiens.net")
    print("graph created!")
    prop = Propagater(ppi, 0.9)
    prop.propagate({1, 2, 100})
    print("propagation complete!")
    jh = 76