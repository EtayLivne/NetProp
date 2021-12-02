import json
from results_models import PropagationResultModel
from sortedcontainers import SortedList
from pathlib import Path
from result_analysis.utils import NORM_FUNCTIONS, propagation_to_node_subset


def sorted_randomized_propagation_scores(randomized_results_dir, liquids_norm="l1"):
    sorted_scores = dict()
    for file in Path(randomized_results_dir).glob("*.json"):
        nodes = PropagationResultModel.parse_file(Path(file)).nodes
        for node, data in nodes.items():
            if node not in sorted_scores:
                sorted_scores[node] = SortedList()
            sorted_scores[node].add(propagation_to_node_subset(nodes, [node], liquids_norm=liquids_norm))

    return sorted_scores


def get_rdpn_p_values(result_file, randomized_results_dir):
    result_under_test = PropagationResultModel.parse_file(result_file)
    sorted_randomized_scores = sorted_randomized_propagation_scores(randomized_results_dir)
    nodes_under_test = result_under_test.nodes
    p_values = dict()
    bad_nodes = 0
    for node in nodes_under_test:
        # if len(sorted_randomized_scores.get(node, [])) < len(sorted_randomized_scores):
        if node not in sorted_randomized_scores:
            print(f"node id {node} does not appear in all randomizations - dropping it")
            bad_nodes += 1
        score = propagation_to_node_subset(nodes_under_test, [node])
        sorted_randomized_scores[node].add(score)
        p_values[node] = 1 - sorted_randomized_scores[node].index(score) / len(sorted_randomized_scores[node])

    print(f"bad nodes: {bad_nodes}")

    return p_values


def record_significant_nodes(result_file, randomized_results_dir, output_path, p_value_threshold=0.05):
    p_values = get_rdpn_p_values(result_file, randomized_results_dir)
    with open(output_path, 'w') as handler:
        json.dump({node: p_value for node, p_value in p_values.items() if p_value < p_value_threshold}, handler,
                  indent=4)


if __name__ == "__main__":
    record_significant_nodes(r"C:\studies\code\NetProp\src\temp\propagations\h_sapiens_from_covid\merged_covid_source.json",
                             r"C:\studies\code\NetProp\src\temp\propagations\h_sapiens_from_covid\merged_covid_source_randomized",
                             r"C:\studies\code\NetProp\src\temp\propagations\h_sapiens_from_covid\significant_nodes.json")