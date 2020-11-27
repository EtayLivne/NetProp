import networkx as nx
from copy import deepcopy
from itertools import chain
from math import sqrt, pow, floor
from collections.abc import Iterable

from datastructures.graph_datastructures import Protein, construct_prior_set
from data_extractors import acquire_prior_set_data
# TODO finish refactor of collection of functions to class
# TODO maybe save both normalized and unnormalized edge weights?
# TODO support multiprocessing, shouldn't be difficult (famous last words, and this is really more of a bonus)
WEIGHT = "weight"
DATA = 2



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
    _DEFUALT_HALT_CONDITION_GAP = 10e-2
    _DEFAULT_LIQUID_TYPE = "__default_liquid"

    def __init__(self, graph, confidence_coefficient, halt_condition_gap=_DEFUALT_HALT_CONDITION_GAP):
        self.graph = graph
        self.confidence_coefficient = confidence_coefficient
        self.halt_condition_gap = halt_condition_gap

        self._normalize_edge_weights()
        self._reset_propagation_parameters()

    @staticmethod
    def validate_propagation_node_sets(network, *nodes):
        missing_nodes = [node.id for node in chain(*nodes) if node not in network.nodes.keys()]
        if missing_nodes:
            print(f'The following nodes are specified in propagation node_set parameters but are not present in '
                  f'the propagation network: {missing_nodes}')

    def _reset_propagation_parameters(self, liquid_types=None):
        for node in self.graph:
            self._reset_liquid(node, liquid_types=liquid_types)

    # assumes graph is directed, otherwise each edge will be normalized twice, which is bad!
    def _normalize_edge_weights(self):
        for node in self.graph:
            edges = [self.graph[node][neighbor] for neighbor in self.graph.neighbors(node)]
            sum_of_weights = sum(edge[WEIGHT] for edge in edges)
            for edge in edges:
                edge[WEIGHT] = float(edge[WEIGHT]) / sum_of_weights

    def __get_node_properties(self, node):
        return self.graph.nodes[node]

    def _get_prev_liquid(self, node, liquid_type):
        node_properties = self.__get_node_properties(node)
        return node_properties[self._PREV_LIQUID_PROPERTY_KEY].get(liquid_type, 0)

    def _set_prev_liquid(self, node, prev_liquid, liquid_type):
        node_properties = self.__get_node_properties(node)
        node_properties[self._PREV_LIQUID_PROPERTY_KEY][liquid_type] = prev_liquid

    def _reset_liquid(self, node, liquid_types=None):
        node_properties = self.__get_node_properties(node)
        if liquid_types:
            for liquid_type in liquid_types:
                node_properties[self._PREV_LIQUID_PROPERTY_KEY][liquid_type] = 0
                node.liquids[liquid_type] = 0
        else:
            node_properties[self._PREV_LIQUID_PROPERTY_KEY] = dict()
            node.liquids = dict()

    # Stores a snapshot of current liquid amount in each node into the node's prev_liquid property
    def _snapshot_to_prev_liquid(self, node, liquid_types):
        node_prev_liquids = self.__get_node_properties(node)[self._PREV_LIQUID_PROPERTY_KEY]
        try:
            node_prev_liquids.update({liquid_type: node.liquids[liquid_type] for liquid_type in liquid_types})
        except:
            jhgf = 0

    def _propagate_once(self, subgraph, prior_set_subgraph, liquid_types):
        alpha = self.confidence_coefficient
        for node in subgraph:
            self._snapshot_to_prev_liquid(node, liquid_types=liquid_types)

        # If liquids dict is empty in any node in prior set then this is the first propagation step

        for node in subgraph.nodes():
            for neighbor in subgraph.neighbors(node):
                edge_capacity = subgraph[node][neighbor][WEIGHT] * alpha
                for liquid_type in liquid_types:
                    neighbor.liquids[liquid_type] = neighbor.liquids.get(liquid_type, 0) +\
                                                    node.liquids.get(liquid_type, 0) * edge_capacity
            for liquid_type in liquid_types:
                node.liquids[liquid_type] = node.liquids.get(liquid_type, 0) - self._get_prev_liquid(node, liquid_type)
            # pump new liquid for next round

        for node in prior_set_subgraph:
            node.liquids.update({liquid_type: node.liquids[liquid_type] + 1 - alpha for liquid_type in node.source_of})


    def propagate(self, prior_set):
        self.validate_propagation_node_sets(self.graph, prior_set)

        liquid_types = set.union(*[prior.source_of for prior in prior_set])
        self._reset_propagation_parameters(liquid_types=liquid_types)

        prior_set_subgraph = nx.subgraph(self.graph, [node for node in self.graph.nodes.keys() if node in prior_set])
        for node in prior_set_subgraph.nodes():
            node.source_of = next(prior for prior in prior_set if prior.id == node.id).source_of

        discovered_subgraph_nodes = set(prior_set_subgraph.nodes())
        newest_subgraph_nodes = set(prior_set_subgraph.nodes())

        # set gaps for all nodes in prior set to >> halt condition
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
    x = acquire_prior_set_data()
    prior_set_data = {id: data.infection_roles for id, data in acquire_prior_set_data().items()}
    ppi = graph_from_file(r"..\data\H_sapiens.net")
    print("graph created!")

    prop = Propagater(ppi, 0.3)
    prop.propagate(construct_prior_set(prior_set_data))
    print("propagation complete!")
    top_k = sorted(prop.graph.nodes(), key=lambda node: sqrt(sum(pow(node.liquids[liquid], 2) for liquid in node.liquids)), reverse=True)
    print("done")