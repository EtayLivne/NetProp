from math import sqrt
from typing import Dict
from results_models import PropagationNodeModel


def l1_norm(nums):
    if not isinstance(nums, list):
        nums = list(nums)
    return sum(nums)


def l2_norm(nums):
    return sqrt(sum([num ** 2 for num in nums]))


NORM_FUNCTIONS = {
    "l1": l1_norm,
    "l2": l2_norm
}


class ProcessedIterator:
    """
    A list iterator that yields the result of activating a function on each element rather than the elements themselves
    """

    def __init__(self, raw_items, processing_func, processing_args, processing_kwargs):
        self._raw_items = raw_items
        self._func = processing_func
        self._args = processing_args
        self._kwargs = processing_kwargs

    def __iter__(self):
        return self

    def __next__(self):
        for item in self._raw_items:
            yield self._func(item, *self._args, **self._kwargs)


def filter_results(nodes: Dict[str, PropagationNodeModel], liquids_filter, nodes_filter):
    for node_name, node in nodes.items():
        if nodes_filter(node_name):
            node.liquids = {liquid: score for liquid, score in node.liquids.items() if liquids_filter(liquid)}
        else:
            nodes.pop(node_name)


def propagation_to_node_subset(nodes: Dict[str, PropagationNodeModel], subset_node_names,
                               nodes_norm="l1", liquids_norm="l1"):
    return NORM_FUNCTIONS[nodes_norm]([NORM_FUNCTIONS[liquids_norm](list(nodes[n].liquids.values())) for
                                       n in subset_node_names])