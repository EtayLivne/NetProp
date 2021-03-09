from math import sqrt
from copy import deepcopy
from itertools import chain

import networkx as nx

from datastructures.graph_datastructures import Protein


class Propagater:
    _PREV_LIQUIDS_PROPERTY_KEY = "prev_liquids"
    _DEFUALT_HALT_CONDITION_GAP = 10e-2
    _DEFAULT_LIQUID_TYPE = "__default_liquid"
    _WEIGHT = 2
    _NX_EDGE_WEIGHT_KEY = "weight"

    def __init__(self, graph, source_condifence_coef, halt_condition_gap=_DEFUALT_HALT_CONDITION_GAP):
        self.graph = graph
        self.source_confidence_coef = source_condifence_coef
        self.halt_condition_gap = halt_condition_gap

        self._normalize_edge_weights()
        self._reset_liquids()

    def gene_knockout_copy(self, gene_knockout_ids):
        knockout_graph = deepcopy(self.graph)
        knockout_graph.remove_nodes_from({Protein(id=gene_id) for gene_id in gene_knockout_ids})
        return Propagater(knockout_graph, self.source_confidence_coef, halt_condition_gap=self.halt_condition_gap)

    def validate_propagation_node_sets(self, *nodes):
        missing_nodes = [node.id for node in chain(*nodes) if node not in self.graph.nodes.keys()]
        if missing_nodes:
            print(f'The following nodes are specified in propagation node_set parameters but are not present in '
                  f'the propagation network: {missing_nodes}')
            return False
        return True

    # assumes graph is directed, otherwise each edge will be normalized twice, which is bad!
    def _normalize_edge_weights(self):
        for node in self.graph:
            edges = [self.graph[node][neighbor] for neighbor in self.graph.neighbors(node)]
            sum_of_weights = sum(edge[self._NX_EDGE_WEIGHT_KEY] for edge in edges)
            for edge in edges:
                edge[self._NX_EDGE_WEIGHT_KEY] = float(edge[self._NX_EDGE_WEIGHT_KEY]) / sum_of_weights

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
                edge_capacity = subgraph[node][neighbor][self._NX_EDGE_WEIGHT_KEY] * alpha
                for liquid_type in liquid_types:
                    neighbor.liquids[liquid_type] = neighbor.liquids.get(liquid_type, 0) + \
                                                    self._get_prev_liquid(node, liquid_type) * edge_capacity
            for liquid_type in liquid_types:
                node.liquids[liquid_type] = node.liquids.get(liquid_type, 0) - self._get_prev_liquid(node, liquid_type)
            # pump new liquid for next round

        for node in prior_set_subgraph:
            node.liquids.update({liquid_type: node.liquids[liquid_type] + 1 - alpha for liquid_type in node.source_of})


    def propagate(self, prior_set, normalize_flow=False):
        if not self.validate_propagation_node_sets(self.graph, prior_set):
            raise ValueError("cannot propagate because prior_set refers to unexisting nodes")

        

        liquid_types = set.union(*[prior.source_of for prior in prior_set])
        self._reset_liquids(liquid_types=liquid_types)

        prior_set_subgraph = nx.subgraph(self.graph, [node for node in self.graph.nodes.keys() if node in prior_set])
        #TODO get rid of this - source_of and target_of are now part of the Protein package and will have been updated seperately.
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

            if len(self.graph) > len(discovered_subgraph):
                newest_subgraph_nodes = set.union(*[set(self.graph.neighbors(node)) for node in newest_subgraph_nodes])
                discovered_subgraph_nodes.update(newest_subgraph_nodes)

        if normalize_flow:
            for liquid_type in liquid_types:
                type_max = max([node.liquids.get(liquid_type, 0) for node in self.graph.nodes if node not in prior_set])
                for node in self.graph.nodes:
                    node.liquids[liquid_type] = float(node.liquids.get(liquid_type, 0)) / type_max