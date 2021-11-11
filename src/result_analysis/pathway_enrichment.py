import scipy.stats as sp_stats
from math import sqrt
import os
import json
from common.network_loaders import HSapeinsNetworkLoader

from results_models import PropagationResultModel

P_VALUE = 1

def l1_norm(nums):
    return sum(nums)

def l2_norm(nums):
    return sqrt(sum([num**2 for num in nums]))

norm_functions = {
    "l1": l1_norm,
    "l2": l2_norm
}

def pathway_enrichment_p_value(pathway, propagation_scores_a, propagation_scores_b):
    paired_tests = [propagation_scores_a[protein]["rank"] - propagation_scores_b[protein]["rank"] for protein in pathway]
    return 1 if all([diff == 0 for diff in paired_tests]) else sp_stats.wilcoxon(paired_tests)[P_VALUE]


def calculate_ranks(propagated_nodes: dict, norm="l1"):
    norm_func = norm_functions.get(norm, None)
    if not norm_func:
        raise ValueError(f"unrecognized norm function {norm}. Available options are: {norm_functions.keys()}")
    ranking = sorted(list(propagated_nodes.keys()),
                     key=lambda node: norm_func(propagated_nodes[node].get("liquids", {}).values()))
    for i in range(len(ranking)):
        propagated_nodes[ranking[i]]["rank"] = i


def pathway_enrichment(pathways: dict, reference_propagation: dict, propagation_to_analyze: dict,
                       reference_already_ranked=True, norm="l1"):
    calculate_ranks(propagation_to_analyze, norm=norm)
    if not reference_already_ranked:
        calculate_ranks(reference_propagation, norm=norm)

    return {pathway: pathway_enrichment_p_value(genes, propagation_to_analyze, reference_propagation) for
            pathway, genes in pathways.items()}


def load_pathways(pathway_file: str):
    """
    structure of a line is pathway__name link_to_webpage pathway_genes
    :param pathway_file:
    :return:
    """
    pathways = dict()
    with open(pathway_file, 'r') as handler:
        for line in handler.readlines():
            as_list = line.split()
            pathways[as_list[0]] = as_list[2:]

    return pathways


def extract_node_dict(propagation_results: PropagationResultModel, relevant_liquids: dict=None):
    nodes = dict()
    relevant_liquids = relevant_liquids or []
    for node in propagation_results.nodes:
        try:
            int(node.id)
        except ValueError:
            continue
        nodes[node.id] = {
            "liquids": {k: v for k, v in node.liquids.items() if k in relevant_liquids or node.liquids}
        }
    return nodes


def keep_only_parsable_pathwaeys(pathways, nodes):
    all_pathway_genes = set.union(*[set(v) for v in pathways.values()])
    to_remove = []
    for pathway, genes in pathways.items():
        if [g for g in genes if g not in nodes]:
            to_remove.append(pathway)

    for pathway in to_remove:
        pathways.pop(pathway)


def global_enrichment_map(parent_dir, threshold_p_value=0.05):
    pathway_to_knockout = {}
    knockout_to_pathway = {}
    for file in os.listdir(parent_dir):
        knockout = file.split("_")[0]
        with open(os.path.join(parent_dir, file), 'r') as handler:
            knockout_enrichment = json.load(handler)
            for pathway, score in knockout_enrichment.items():
                is score <= threshold_p_value:
                pathway_to_knockout



def main():
    pathways = load_pathways(r"C:\studies\code\NetProp\data\c2.cp.v7.4.entrez.gmt")
    reference_propagation_path = r"C:\studies\code\NetProp\src\temp\propagations\h_sapiens_from_covid\all_covid_sources.json"
    reference_propagation_nodes = extract_node_dict(PropagationResultModel.parse_file(reference_propagation_path))
    calculate_ranks(reference_propagation_nodes)


    propagations_dir_path = r'C:\studies\code\NetProp\src\temp\propagations\gene_knockouts'
    for propagation_result_file in os.listdir(propagations_dir_path):
        prop_result = PropagationResultModel.parse_file(os.path.join(propagations_dir_path, propagation_result_file))
        extracted_nodes = extract_node_dict(prop_result)
        keep_only_parsable_pathwaeys(pathways, extracted_nodes)
        pathway_enrichment_map = pathway_enrichment(pathways, reference_propagation_nodes, extracted_nodes)
        with open(os.path.join(r'C:\studies\code\NetProp\src\temp\knockout_pathways_enrichment',
                               propagation_result_file),
                  "w") as handler:
            json.dump(pathway_enrichment_map, handler, indent=4)
        # os.remove(os.path.join(propagations_dir_path, propagation_result_file))


if __name__ == "__main__":
    main()