from itertools import chain
from contextlib import contextmanager
from functools import reduce, cached_property

import numpy as np
import scipy as sp
import pandas as pd
import scipy.linalg
import networkx as nx
import scipy.sparse as sparse
from scipy.sparse.linalg import inv
import numpy.lib.scimath as np_scimath


class Propagater:
    # _LIQUIDS_ATTR = "liquids"
    _NODES_COLUMN = "nodes"

    def __init__(self, network: nx.Graph, source_confidence: float,
                 prior_set: list[str]=None, suppressed_nodes: list[str]=None,
                 min_gap: float=None, max_steps: int=None):
        self._network = None
        self._prior_set = None
        self._suppressed_nodes = None
        self._source_confidence = None
        self._min_gap = None
        self._max_steps = None

        self.network = network
        self.prior_set = prior_set
        self.suppressed_nodes = suppressed_nodes or []
        self.source_confidence = source_confidence
        self.min_gap = min_gap
        self.max_steps = max_steps

        self._state = pd.DataFrame({self._NODES_COLUMN: list(self.network.nodes())})
        self._state = self._state.set_index(self._NODES_COLUMN)
        self._unsuppressed_network = None

    @property
    def suppressed_nodes(self):
        return self._suppressed_nodes

    @suppressed_nodes.setter
    def suppressed_nodes(self, suppressed_nodes: list[str]):
        self._validate_node_set(suppressed_nodes)
        self._suppressed_nodes = suppressed_nodes

    @property
    def prior_set(self):
        return self._prior_set

    @prior_set.setter
    def prior_set(self, prior_set):
        self._validate_node_set(prior_set)
        self._prior_set = prior_set

    @property
    def network(self):
        return self._network

    @network.setter
    def network(self, network):
        if not isinstance(network, nx.Graph):
            raise TypeError(f"May only set the network attribute of a propagater to an object of type nx.Graph,"
                            f"got object of type {type(network)} instead")

        self._network = network
        self._reset_weight_matrix()

    @property
    def min_gap(self):
        return self._min_gap

    @min_gap.setter
    def min_gap(self, min_gap):

        try:
            if min_gap is not None:
                min_gap = float(min_gap)
                if min_gap < 0:
                    raise ValueError
        except ValueError:
            raise TypeError(f"min_gap must be a non negative number, received {type(min_gap)} with value {min_gap}")

        self._min_gap = min_gap

    @property
    def max_steps(self):
        return self._max_steps

    @max_steps.setter
    def max_steps(self, max_steps):
        try:
            if max_steps is not None:
                max_steps = int(max_steps)
                if max_steps < 1:
                    raise ValueError
        except ValueError:
            raise TypeError(f"max_steps must be an integer, received {type(max_steps)} with value {max_steps}")

        self._max_steps = max_steps

    def propagate(self, prior_set: list[str]=None, suppressed_nodes: list[str]=None, source_confidence: float=None,
                  method="iterative"):
        self._validate_and_update_params(prior_set, suppressed_nodes, source_confidence)
        # with self._suppress_nodes():  # The day will come... I think
        suppressed_indexes = self._state[list(self.suppressed_nodes)] if self.suppressed_nodes else []
        for liquid, liquid_prior_vector in self.get_liquid_specific_vectors():
            liquid_prior_vector = self.source_confidence * liquid_prior_vector
            if method == "iterative":
                state_vector = self._iterative_propagation(liquid_prior_vector, suppressed_indexes)
            else:
                state_vector = self._analytical_propagation(liquid_prior_vector)
            self._state[liquid] = np.array(state_vector).flatten()

    def get_network_state(self, by_liquid: str=None):
        if by_liquid:
            try:
                return self._state[by_liquid]
            except KeyError:
                raise ValueError(f"Propagater has no data on liquid {by_liquid}")
        return self._state

    def get_node_state(self, node):
        return self._state.loc[node]

    def _iterative_propagation(self, prior_vector: sparse.csr_matrix, suppressed_indexes: list[int]):
        state_vector = prior_vector.todense()
        step_counter = self.max_steps
        while True:
            prev_state_vector = state_vector.copy()
            # Note that both weight matrix and prior vector are weighted by prior confidence before this method
            state_vector = self._weight_matrix.dot(prev_state_vector) + prior_vector
            for i in suppressed_indexes:
                state_vector[i, 0] = 0
            if self.min_gap is not None and sp.linalg.norm(state_vector - prev_state_vector) < self.min_gap:
                return state_vector
            if self.max_steps is not None:
                step_counter -= 1
                if step_counter <= 0:
                    break
        return state_vector

    def _analytical_propagation(self, prior_vector: sparse.csr_matrix):
        I = sp.sparse.diags([1] * len(self.network.nodes), format="csr")
        return sp.sparse.linalg.inv(I - self._weight_matrix).dot(prior_vector)

    def _validate_and_update_params(self, prior_set, suppressed_nodes, source_confidence):
        self.prior_set = prior_set or self.prior_set
        self.suppressed_nodes = suppressed_nodes or self.suppressed_nodes
        self.source_confidence = source_confidence or self.source_confidence
        missing = {
            "prior_set": (self.prior_set is None) or len(self.prior_set) == 0,
            "suppressed_nodes": self.suppressed_nodes is None,
            "source_confience": self.source_confidence is  None,
            "network": self.network is None
        }

        if missing_params := [p for p in missing if missing[p]]:
            raise ValueError(f"Failed to initiate propagation because"
                             f" the following params are not set: {missing_params}")

        if len(self.prior_set) == 0:
            raise ValueError("Failed to initiate propagation because prior_set is empty")

        if self.max_steps is None and self.min_gap is None:
            raise ValueError("No halt condition (max step or min gap) set for propagation")

    def get_liquid_specific_vectors(self):
        # Automatically infers which liquids to propagate based on "source of" attr of nodes in prior set
        # Returns: tuples of (liquid_name, p0 vector for prop with that prior set
        liquid_specific_prior_vectors = {}
        for prior in self.prior_set:
            prior_index = self._state.index.get_loc(prior)
            for liquid in self.network.nodes[prior].get("source_of", set()):
                if liquid not in liquid_specific_prior_vectors:
                    liquid_specific_prior_vectors[liquid] = []
                liquid_specific_prior_vectors[liquid].append(prior_index)

        # Transform lists of indexes to sparse column vector of unit indicators for each liquid
        to_return = []
        for liquid, row_indexes in liquid_specific_prior_vectors.items():
            data = np.ones((len(row_indexes)))
            col_indexes = np.zeros((len(row_indexes)))
            to_return.append(
                (
                 liquid,
                 sparse.csr_matrix((data, (row_indexes, col_indexes)),
                                   shape=(len(self._state), 1))
                 )
            )
        return to_return

    @contextmanager
    def _suppress_nodes_but_not_now_not_really_right(self):
        self._unsuppressed_network = self.network
        self.network = nx.subgraph(self.network, [n for n in self.network.nodes if n not in self.suppressed_nodes])

        try:
            yield
        finally:
            self.network = self._unsuppressed_network

    @cached_property
    def _weight_matrix(self) -> sparse.csr_matrix:
        adjacency_matrix = nx.to_scipy_sparse_matrix(self.network)
        weighted_deg_matrix = sp.sparse.diags(1 / np_scimath.sqrt(adjacency_matrix.sum(0).A1), format="csr")
        confidence = 1 - self.source_confidence
        weight_matrix = confidence * weighted_deg_matrix * adjacency_matrix * weighted_deg_matrix
        return weight_matrix

    def _reset_weight_matrix(self):
        if hasattr(self, "_weight_matrix"):
            delattr(self, "_weight_matrix")

    def _validate_node_set(self, node_set):
        try:
            copied_list = list(node_set)
        except TypeError as te:
            raise TypeError(f"suppressed_nodes must be an iterable, got input of type {type(node_set)}") from te

        if unexsisting_nodes := [n for n in copied_list if n not in self.network]:
            raise ValueError(f"Failed to set an input list as suppressed nodes for propagation because the following"
                             f"nodes aren't in the network: {unexsisting_nodes}")