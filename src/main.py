import networkx as nx
from copy import deepcopy
from math import sqrt, pow

# TODO finish refactor of collection of functions to class
# TODO maybe save both normalized and unnormalized edge weights?
# TODO support multiprocessing, shouldn't be difficult (famous last words, and this is really more of a bonus)
WEIGHT = "weight"


def graph_from_file(file_path):
    source_node_index = 0
    target_node_index = 1
    edge_weight_index = 2

    graph = nx.DiGraph()
    with open(file_path, 'r') as handler:
        for line in handler.readlines():
            values = line.split()
            graph.add_edge(values[source_node_index], values[target_node_index], weight=values[edge_weight_index])


def get_propagater(source, from_file=True):
    graph = graph_from_file(source) if from_file else deepcopy(source)
    return Propagater(graph)


class Propagater:
    def __init__(self, graph):
        self.graph = graph
        self._liquid_snapshot = None
        self.gap = None

        self._normalize_edge_weights()
        self._reset_propagation_parameters()

    def _reset_propagation_parameters(self):
        self._liquid_snapshot = dict()
        self.gap = 0
        for node in self.graph.nodes():
            node.liquid = 0

    # assumes graph is directed, otherwise each edge will be normalized twice, which is bad!
    def _normalize_edge_weights(self):
        for node in self.graph:
            sum_of_weights = sum(edge[WEIGHT] for edge in self.graph.out_edges(node))
            for edge in self.graph.out_edges(node):
                edge[WEIGHT] = float(edge[WEIGHT]) / sum_of_weights

    def _propagate_once(self, subgraph, confidence_coefficient, prior_set):
        alfa = confidence_coefficient
        self._liquid_snapshot = {node: node.liquid for node in subgraph.nodes()}

        if next(iter(prior_set)).liquid == 0:
            for node in prior_set:
                node.liquid = 1 - alfa
                self.gap = sqrt(pow((1-alfa), 2) * len(prior_set))
        else:
            for node in subgraph.nodes():
                if node in prior_set:
                    node.liquid = node.liquid + (1 - alfa)
                for neighbour in subgraph.neighbours(node):
                    neighbour.liquid = neighbour.liquid + node.liquid * subgraph[node][neighbour][WEIGHT] * alfa
                node.liquid = node.liquid - self._liquid_snapshot[node]
            self.gap = sqrt(sum(pow(node.liquid - self._liquid_snapshot[node], 2) for node in subgraph.nodes()))

    def propagate(self, confidence_coefficient, prior_set, halt_condition_gap=10e-4):
        self._reset_propagation_parameters()
        discovered_subgraph_nodes = prior_set
        newest_subgraph_nodes = prior_set

        self.gap = 10*halt_condition_gap  # only important detail is that gap > halt_condition_gap
        while self.gap > halt_condition_gap:
            discovered_subgraph = nx.subgraph(self.graph, discovered_subgraph_nodes)
            self._propagate_once(discovered_subgraph, confidence_coefficient, prior_set)
            newest_subgraph_nodes = set.union(*{set(self.graph.neighbours(node)) for node in newest_subgraph_nodes})
            discovered_subgraph_nodes.union(newest_subgraph_nodes)

    def top_k_candidates(self, k):
        return sorted(self.graph.nodes(), reverse=True)[:k]
