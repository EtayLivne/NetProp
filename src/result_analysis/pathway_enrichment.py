import scipy.stats as sp_stats
from math import sqrt
import os
import json
from common.network_loaders import HSapeinsNetworkLoader
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

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





def global_enrichment_map(parent_dir, threshold_p_value=0.05):
    pathway_to_knockout = {}
    knockout_to_pathway = {}
    for file in os.listdir(parent_dir):
        knockout = file.split("_")[0]
        with open(os.path.join(parent_dir, file), 'r') as handler:
            knockout_enrichment = json.load(handler)
            for pathway, score in knockout_enrichment.items():
                if score <= threshold_p_value:
                    if pathway not in pathway_to_knockout:
                        pathway_to_knockout[pathway] = []
                    pathway_to_knockout[pathway].append(knockout)

                    if knockout not in knockout_to_pathway:
                        knockout_to_pathway[knockout] = []
                    knockout_to_pathway[knockout].append(pathway)

    return pathway_to_knockout, knockout_to_pathway


def ordered_json_dump(dct, key_sorting_function, output_path, reverse_sort=False):
    sorted_keys = sorted(list(dct.keys()), key=key_sorting_function, reverse=reverse_sort)
    with open(output_path, 'w') as handler:
        json.dump({key: dct[key] for key in sorted_keys}, handler, indent=4)


def map_enrichment(data_dir, output_dir, threshold_p_value=0.05):
    pathway_to_knockout, knockout_to_pathway = global_enrichment_map(data_dir, threshold_p_value=threshold_p_value)

    print(f"knockouts that enrich: {len(knockout_to_pathway)}. \n",
          f"average number of pathways enriched: {sum([len(v) for v in knockout_to_pathway.values()]) / len(knockout_to_pathway)}.\n"
          f"Pathways that are enriched: {len(pathway_to_knockout)}.\n",
          f"average enrichment sources: {sum([len(v) for v in pathway_to_knockout.values()]) / len(pathway_to_knockout)}."
          )
    # ordered_json_dump(knockout_to_pathway, lambda knockout: len(knockout_to_pathway[knockout]),
    #                   os.path.join(output_dir, "knockout_to_pathway.json"), reverse_sort=True)
    #
    # ordered_json_dump(pathway_to_knockout, lambda pathway: len(pathway_to_knockout[pathway]),
    #                   os.path.join(output_dir, "pathway_to_knockout.json"), reverse_sort=True)





def main():
    # map_enrichment(r"C:\studies\thesis\code\NetProp\src\temp\knockout_pathways_enrichment",
    #                r'C:\studies\thesis\code\NetProp\src\temp', threshold_p_value=0.00000794912)
    pathways = load_pathways(r"C:\studies\thesis\code\NetProp\data\canonical_pathways.gmt")
    num_bins = 100
    n, bins, patches = plt.hist([len(v) for v in pathways.values()], num_bins, facecolor='blue', alpha=0.5)
    plt.show()
    # pathways = load_pathways(r"C:\studies\code\NetProp\data\c2.cp.v7.4.entrez.gmt")
    # reference_propagation_path = r"C:\studies\code\NetProp\src\temp\propagations\h_sapiens_from_covid\all_covid_sources.json"
    # reference_propagation_nodes = extract_node_dict(PropagationResultModel.parse_file(reference_propagation_path))
    # calculate_ranks(reference_propagation_nodes)
    #
    #
    # propagations_dir_path = r'C:\studies\code\NetProp\src\temp\propagations\gene_knockouts'
    # for propagation_result_file in os.listdir(propagations_dir_path):
    #     prop_result = PropagationResultModel.parse_file(os.path.join(propagations_dir_path, propagation_result_file))
    #     extracted_nodes = extract_node_dict(prop_result)
    #     keep_only_parsable_pathwaeys(pathways, extracted_nodes)
    #     pathway_enrichment_map = pathway_analysis(pathways, reference_propagation_nodes, extracted_nodes)
    #     with open(os.path.join(r'C:\studies\code\NetProp\src\temp\knockout_pathways_enrichment',
    #                            propagation_result_file),
    #               "w") as handler:
    #         json.dump(pathway_enrichment_map, handler, indent=4)
    #     # os.remove(os.path.join(propagations_dir_path, propagation_result_file))


if __name__ == "__main__":
    main()