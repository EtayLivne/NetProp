from typing import Set
from itertools import chain

import numpy as np
import scipy as sp
import scipy.linalg
import networkx as nx
# from pydantic import Field
from scipy.sparse.linalg import inv
import numpy.lib.scimath as np_scimath
# from pydantic.dataclasses import dataclass

from netprop.networks.propagation_network import PropagationContainer, PropagationNetwork
# @dataclass
# class PropagationContainer:
#     source_of: Set[str] = Field(default_factory=set)
#     # target_of: Set[str] = Field(default_factory=set)  #TODO delete if we officially give up on target sets


# class PropagationNetwork(nx.Graph):
#     CONTAINER_KEY = "container"
#     EDGE_WEIGHT = "weight"
#
#     @classmethod
#     def node_has_propagation_container(cls, node_properties):
#         return isinstance(node_properties.get(cls.CONTAINER_KEY, None), PropagationContainer)
#
#     def is_ready_for_propagation(self):
#         nodes_condition = all(self.node_has_propagation_container(data) for node, data in self.nodes(data=True))
#         if not nodes_condition:
#             return False
#         edges_condition = all(data.get("weight", -1) > 0 for _, _, data in self.edges(data=True))
#
#         return nodes_condition and edges_condition


class PropagationError(BaseException):
    pass

class PropagationNetworkError(BaseException):
    pass

class Propagater:
    _PREV_LIQUIDS = "prev_liquids"
    _LIQUIDS = "liquids"
    _EDGE_WEIGHT = PropagationNetwork.EDGE_WEIGHT
    _DEFAULT_HALT_CONDITION_GAP = 10e-5
    NO_MAX_STEPS = -1
    NO_MIN_GAP = -1
    ITERATIVE = "iterative"
    ANALYTICAL = "analytical"


    @staticmethod
    def _edge_degree_coef(node1_deg, node2_deg):
        return float(2) / (node1_deg + node2_deg)

    @classmethod
    def propagation_related_node_metadata(cls):
        return [cls._LIQUIDS, cls._PREV_LIQUIDS, PropagationNetwork.CONTAINER_KEY]

    def __init__(self, network, source_confidence, min_gap=NO_MIN_GAP, max_steps=NO_MAX_STEPS):
        self.network = network
        self.source_confidence = source_confidence
        self.max_steps = max_steps
        self.min_gap = min_gap


        self._reset_liquids()
        self._edges_are_normalized = False
        self._unsuppressed_network = None
        self._propagation_methods = {self.ITERATIVE: self._iterative_propagation,
                                     self.ANALYTICAL: self._analytic_propagation}
        self._adjacency_matrix = self._normalize_edge_weights(nx.to_scipy_sparse_matrix(self.network))

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

    def propagate(self, prior_set, suppressed_set=None, propagation_method="iterative"):
        # validate input
        if not self.validate_propagation_node_sets(self._network, prior_set):
            raise ValueError("Cannot propagate because prior_set refers to nodes that do not exist")

        if self._max_steps == self.NO_MAX_STEPS and self._min_gap == self.NO_MIN_GAP:
            raise PropagationError("No halt condition specified")

        if propagation_method not in self._propagation_methods:
            raise ValueError(f"propagation method {propagation_method} unrecognized. "
                             f"Choose one of the following: {self._PROPAGATION_METHODS}")

        if propagation_method == self.ANALYTICAL and len(suppressed_set) > 0:
            raise ValueError("analytical method with suppressed nodes is not supported in this version :(")

        # Prepare propagation based on parameters
        # self._suppress_nodes(suppressed_nodes=suppressed_set) # new suppression method leaves the nodes in the graph
        # adjacency_matrix = nx.to_scipy_sparse_matrix(self.network)    # since the matrix is now the same in all iterations, do this once at the object level
        # adjacency_matrix = self._normalize_edge_weights(self._adjacency_matrix) # since the matrix is now the same in all iterations, do this once at the object level
        liquids = {}
        for prior in prior_set:
            for liquid in self.network.nodes[prior][self.network.CONTAINER_KEY].source_of:
                if liquid in liquids:
                    liquids[liquid].add(prior)
                else:
                    liquids[liquid] = {prior}
        self._reset_liquids(liquid_types=liquids.keys())

        # Propagate
        for liquid, priors in liquids.items():
            prior_vector = self.source_confidence * np.array([1 if n in priors else 0 for n in self.network.nodes])
            suppressed_indexes = np.nonzero( np.array([1 if n in suppressed_set else 0 for n in self.network.nodes]))[0]
            state_vector = self._propagation_methods[propagation_method](prior_vector, suppressed_indexes)
            for index, node in enumerate(self.network.nodes):
                self.network.nodes[node][self._LIQUIDS][liquid] = state_vector[index]

        # self._suppress_nodes(reverse=True)

    def _iterative_propagation(self, prior_vector, suppressed_indexes):
        state_vector = np.array([])
        step_counter = self.max_steps
        while True:
            prev_state_vector = np.copy(state_vector) if state_vector.any() else np.copy(prior_vector)
            state_vector = self._adjacency_matrix.dot(prev_state_vector) + prior_vector
            for i in suppressed_indexes:
                state_vector[i] = 0

            if self.min_gap != self.NO_MIN_GAP and sp.linalg.norm(state_vector - prev_state_vector) < self.min_gap:
                break
            if self.max_steps != self.NO_MAX_STEPS:
                step_counter -= 1
                if not step_counter:
                    break
        return state_vector

    def _analytic_propagation(self, prior_vector, suppressed_set):
        I = sp.sparse.diags([1] * len(self.network.nodes), format="csr")
        return sp.sparse.linalg.inv(I - self._adjacency_matrix).dot(prior_vector)

    #TODO support multiple forms of normalization (via parameter in propagation config)
    # Only works for symmetric (i.e undirected) networks!
    def _normalize_edge_weights(self, adjacency_matrix, apply_confidence_coef=True):
        weighted_deg_matrix = sp.sparse.diags(1 / np_scimath.sqrt(adjacency_matrix.sum(0).A1), format="csr")
        coef = 1 - self.source_confidence if apply_confidence_coef else 1
        return coef * weighted_deg_matrix * adjacency_matrix * weighted_deg_matrix

    def node_liquids(self, node_id):
        return self.network.nodes[node_id][self._LIQUIDS]

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

    def _suppress_nodes(self, suppressed_nodes=None, reverse=False):

        if reverse and not bool(self._unsuppressed_network):
            raise PropagationError("No nodes to un-suppress")
        if bool(self._unsuppressed_network) and not reverse:
            raise PropagationError("May not suppress nodes while another set of nodes is already supressed")

        if reverse:
            self._network = self._unsuppressed_network
            self._unsuppressed_network = None

        else:
            # nodes that need to be removed are those explicitly specified and those orphaned by their removal.
            suppressed_nodes = [node for node in suppressed_nodes] if suppressed_nodes else list()
            orpahned_neighbors = []
            for node in suppressed_nodes:
                orpahned_neighbors.extend([neighbor for neighbor in self._network.neighbors(node) if
                                           self._network.degree(neighbor) == 1])
            suppressed_nodes.extend(orpahned_neighbors)
            suppressed_network_nodes = [node for node in self._network.nodes if node not in suppressed_nodes]
            # print(f"there are {len(suppressed_network_nodes)} nodes in the subgraph when {suppressed_nodes} are suppressed")

            self._unsuppressed_network = self._network
            try:
                self._network = nx.subgraph(self._network, suppressed_network_nodes)
            except TypeError:
                self._unsuppressed_network = None
                raise TypeError(f"failed to suppress {suppressed_nodes} because it is not an iterable of node ids")




#TODO remove propagation test code below. It's dirty!
# network = PropagationNetwork()
# network.add_edge(1, 2, weight=0.5)
# network.nodes[1][PropagationNetwork.CONTAINER_KEY] = PropagationContainer(source_of={"liquid"})
# network.nodes[2][PropagationNetwork.CONTAINER_KEY] = PropagationContainer(target_of={"liquid"})
# p = Propagater(network, 0.5, max_steps=4)
# p.propagate({1})
# x = 5