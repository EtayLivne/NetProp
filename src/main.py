import networkx as nx
from copy import deepcopy
from math import sqrt, pow, floor

# TODO finish refactor of collection of functions to class
# TODO maybe save both normalized and unnormalized edge weights?
# TODO support multiprocessing, shouldn't be difficult (famous last words, and this is really more of a bonus)
WEIGHT = "weight"
DATA = 2
from data_structures.graph_datastructures import Protein


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
    _PREV_LIQUID_PROPERTY_KEY = "prev_liquid"

    def __init__(self, graph):
        self.graph = graph
        self.gap = None

        self._normalize_edge_weights()
        self._reset_propagation_parameters()

    def _reset_propagation_parameters(self):
        self.gap = 0
        for node in self.graph:
            node.liquid = 0
            self._set_prev_liquid(node, reset=True)

    # assumes graph is directed, otherwise each edge will be normalized twice, which is bad!
    def _normalize_edge_weights(self):
        for node in self.graph:
            edges = [self.graph[node][neighbor] for neighbor in self.graph.neighbors(node)]
            sum_of_weights = sum(edge[WEIGHT] for edge in edges)
            for edge in edges:
                edge[WEIGHT] = float(edge[WEIGHT]) / sum_of_weights

    def __get_node_properties(self, node):
        return self.grap.nodes[node]

    def _get_prev_liquid(self, node, source_node):
        node_properties = self.__get_node_properties(self, node)
        return node_properties[self._PREV_LIQUID_PROPERTY_KEY][source_node]

    def _set_prev_liquid(self, node, prev_liquid, source_node):
        node_properties = self.__get_node_properties(self, node)
        node_properties[self._PREV_LIQUID_PROPERTY_KEY][source_node] = prev_liquid

    def _reset_liquid(self, node, source_nodes=None):
        node_properties = self.__get_node_properties(node)
        if source_nodes:
            for source_node in source_nodes:
                del(node_properties[self._PREV_LIQUID_PROPERTY_KEY][source_node])
        else:
            node_properties[self._PREV_LIQUID_PROPERTY_KEY] = dict()

    def _snapshot_to_prev_liquid(self, node):
        node_prev_liquids = self.__get_node_properties(node)[self._PREV_LIQUID_PROPERTY_KEY]
        node_prev_liquids.update(node.liquids)

    def _propagate_once(self, subgraph, prior_set_subgraph, confidence_coefficient):
        alfa = confidence_coefficient
        for node in subgraph:
            self._snapshot_to_prev_liquid(node)

        if next(iter(prior_set_subgraph)).liquid == 0:
            for node in prior_set_subgraph:
                node.liquid = 1 - alfa
            self.gap = sqrt(pow((1-alfa), 2) * len(prior_set_subgraph))
        else:
            for node in subgraph.nodes():
                for neighbor in subgraph.neighbors(node):
                    neighbor.liquid = neighbor.liquid + self._get_prev_liquid(node) * subgraph[node][neighbor][WEIGHT] * alfa
                node.liquid = node.liquid - self._get_prev_liquid(node)
                # pump new liquid for next round
                if node in prior_set_subgraph:
                    node.liquid = node.liquid + (1 - alfa)

            self.gap = sqrt(sum(pow(node.liquid - self._get_prev_liquid(node), 2) for node in subgraph.nodes()))

    def propagate(self, confidence_coefficient, prior_set_ids, halt_condition_gap=10e-4):
        self._reset_propagation_parameters()
        prior_set = {node for node in self.graph.nodes if node.id in prior_set_ids}
        prior_set_subgraph = nx.subgraph(self.graph, prior_set)

        discovered_subgraph_nodes = set(prior_set)
        newest_subgraph_nodes = set(prior_set)
        self.gap = 10*halt_condition_gap  # only important detail is that gap > halt_condition_gap
        while self.gap > halt_condition_gap:
            discovered_subgraph = nx.subgraph(self.graph, discovered_subgraph_nodes)
            self._propagate_once(discovered_subgraph, prior_set_subgraph, confidence_coefficient)
            newest_subgraph_nodes = set.union(*[set(self.graph.neighbors(node)) for node in newest_subgraph_nodes])
            discovered_subgraph_nodes.update(newest_subgraph_nodes)

    def top_k_candidates(self, k):
        top_k = floor(len(self.graph) / k)
        return sorted(self.graph.nodes(), reverse=True)[:top_k]

if __name__ == "__main__":
    ppi = graph_from_file(r"C:\studies\thesis\code\NetProp\data\H_sapiens.net")
    print("graph created!")
    prop = Propagater(ppi)
    prop.propagate(0.9, {1}, halt_condition_gap=0.01)
    print("propagation complete!")
    top_10 = prop.top_k_candidates(10)
    jh = 76