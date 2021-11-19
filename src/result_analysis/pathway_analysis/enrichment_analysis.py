from .pathways import Pathway, PathwayManager, tag_unknown_genes

import os
from math import sqrt
from typing import List, Set
from sortedcontainers import SortedList
from results_models import PropagationResultModel
from functools import partial


def l1_norm(nums):
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
        self.__func = processing_func
        self._args = processing_args
        self._kwargs = processing_kwargs

    def __iter__(self):
        return self

    def __next__(self):
        for item in self._raw_items:
            yield self._func(item, *self._args, **self._kwargs)


def calc_pathway_propagation(nodes: dict, pathway: Pathway, norm="l1"):
    return NORM_FUNCTIONS[norm]([nodes[n].get("liquids", dict()).values() for n in pathway.genes])


def calc_pathway_p_value(pathway: Pathway,
                         tested_propagation_nodes: dict, randomized_propagation_nodes: List[dict],
                         norm="l1"):
    # initializing the SortedList item from a custom iterator to avoid copying the original list during SortedList construction
    pathway_propagation_score_iterator = ProcessedIterator(randomized_propagation_nodes, [pathway], {"norm": norm})
    sorted_pathway_propagation_scores = SortedList(pathway_propagation_score_iterator)
    test_propagation = calc_pathway_propagation(tested_propagation_nodes, pathway, norm=norm)
    sorted_pathway_propagation_scores.add(test_propagation)
    return 1 - sorted_pathway_propagation_scores.index(test_propagation) / len(sorted_pathway_propagation_scores)


def result_filter(liquids_filter, nodes_filter, propagation_results: PropagationResultModel):
    for node in list(propagation_results.nodes.keys()):
        if nodes_filter(node):
            node.liquids = {liquid: score for liquid, score in node.liquids.items() if liquids_filter(liquid)}
        else:
            node.liquids.pop(node)


def p_value_worker_main(queue,
                        propagation_under_test_file: str,
                        randomized_propagations_dir: str,
                        pathways_file: str,
                        relevant_liquids: Set = None,
                        norm="l1"):
    pathways_manager = PathwayManager.from_file_path(pathways_file)
    prop_under_test = PropagationResultModel.parse_file(propagation_under_test_file)
    nodes = prop_under_test.nodes
    if relevant_liquids:
        liquid_filter = lambda liquid: liquid in relevant_liquids
    else:
        liquid_filter = lambda liquid: True
    node_filter = lambda node: node in nodes
    filter_results = partial(result_filter, liquid_filter, node_filter)
    filter_results(nodes)
    tag_unknown_genes(pathways_manager, prop_under_test.nodes)

    while True:
        request = queue.get(block=True)

        if request == "HALT":
            break
        output_path = randomization_request["output_path"]
        switch_factor = randomization_request["edge_switch_factor"]
        print(f"now handling request for {output_path}")
        network = network_loader.load_network()
        generate_rank_equivalent_random_network(network, switch_factor)
        network_loader_class.record_network(network, output_path)


# assumption: genes in all randomized propagations are a subset of genes in tested propagation (for tag purposes)
def record_pathway_p_values(tested_propagation_result_file: str,
                            randomized_propagations_result_dir: str,
                            pathways_file: str,
                            relevant_liquids: Set = None,
                            norm: str = "l1"):
    pathway_p_values = dict()
    for pathway in pathways_manager.pathways(blacklist=["unknown_genes"]):
        sorted_propagation_values = SortedList([])
        for file in os.listdir(randomized_propagations_result_dir):
            prop = PropagationResultModel.parse_file(os.path.join(randomized_propagations_result_dir, file))
            sorted_propagation_values.add(calc_pathway_propagation(filter_result(prop), pathway, norm=norm))

        tested_prop_val = calc_pathway_propagation(tested_prop.nodes, pathways_manager, norm=norm)
        sorted_propagation_values.add(tested_prop_val)

        pathway_p_values[pathway] = \
            1 - sorted_propagation_values.index(tested_prop.nodes) / len(sorted_propagation_values)

    return pathway_p_values
