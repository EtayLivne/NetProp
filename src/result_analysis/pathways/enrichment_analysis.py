from result_analysis.pathways.pathways import Pathway, PathwayManager, tag_unknown_genes

import json
import os
from math import sqrt
from typing import List, Set, Dict
from sortedcontainers import SortedList
from results_models import PropagationResultModel, PropagationNodeModel
from functools import partial
from multiprocessing import Pool, Process
from pathlib import Path

from result_analysis.utils import filter_results, propagation_to_node_subset


def propagation_to_pathway(pathways: Dict[str, Pathway], norm: str, propagation_result_files: Dict[str, str]):
    pathway_dict = {pathway_name: dict() for pathway_name in pathways}
    for file_name, file_path in propagation_result_files.items():
        nodes = PropagationResultModel.parse_file(file_path).nodes
        for name, pathway in pathways.items():
            pathway_dict[name][file_name] = propagation_to_node_subset(nodes, pathway.genes,
                                                                       nodes_norm=norm, liquids_norm=norm)

    # TODO actual mechanism for temporary files
    with open(os.path.join(r"C:\studies\code\NetProp\src\temporary_dir",
                           f"temp_pathways_{os.getpid()}.json"), 'w') as handler:
        json.dump(pathway_dict, handler)


def load_gene_filtered_pathways(propagation_nodes, pathways_file, relevant_liquids):
    pathways_manager = PathwayManager.from_file_path(pathways_file)
    if relevant_liquids:
        liquid_filter = lambda liquid: liquid in relevant_liquids
    else:
        liquid_filter = lambda liquid: True
    node_filter = lambda node: node in propagation_nodes
    filter_results(propagation_nodes, liquid_filter, node_filter)
    tag_unknown_genes(pathways_manager, propagation_nodes)
    return {name: pathway for name, pathway in pathways_manager.pathways(blacklist={"unknown_genes"})}


# assumption: genes in all randomized propagations are a subset of genes in tested propagation (for tag purposes)
def get_pathway_propagation_scores(propagation_under_test_file: str,
                                   randomized_propagations_result_dir: str,
                                   pathways_file: str,
                                   relevant_liquids: Set = None,
                                   norm: str = "l1"):

    # Prepare pathways
    # pathways_manager = PathwayManager.from_file_path(pathways_file)
    # prop_under_test = PropagationResultModel.parse_file(propagation_under_test_file)
    # nodes = prop_under_test.nodes
    # if relevant_liquids:
    #     liquid_filter = lambda liquid: liquid in relevant_liquids
    # else:
    #     liquid_filter = lambda liquid: True
    # node_filter = lambda node: node in nodes
    # filter_results(nodes, liquid_filter, node_filter)
    # tag_unknown_genes(pathways_manager, prop_under_test.nodes)
    # filtered_pathways = {name: pathway for name, pathway in pathways_manager.pathways(blacklist={"unknown_genes"})}
    prop_under_test = PropagationResultModel.parse_file(propagation_under_test_file)
    filtered_pathways = load_gene_filtered_pathways(prop_under_test.nodes, pathways_file, relevant_liquids)
    # Prepare map of files to iterate
    files_dict = {Path(file).stem.split("_")[-1]: os.path.join(randomized_propagations_result_dir, file) for
                  file in os.listdir(randomized_propagations_result_dir)}
    files_dict["original"] = propagation_under_test_file
    files_dict_keys = list(files_dict.keys())

    # Launch multiprocessed pathway propagation comparison
    processes = []
    for i in range(10):
        partial_files_dict = {k: files_dict[k] for
                              k in files_dict_keys[i*len(files_dict)//10:(i+1)*len(files_dict)//10]}
        p = Process(target=propagation_to_pathway,
                    args=(filtered_pathways, norm, partial_files_dict))
        p.start()
        processes.append(p)

    for p in processes:
        p.join()

    # Create unified dict
    pathway_scores = dict()
    for file in Path(r'C:\studies\code\NetProp\src\temporary_dir').glob("*.json"):
        with open(file, 'r') as handler:
            file_pathways = json.load(handler)
        if not pathway_scores:
            pathway_scores = file_pathways
        else:
            for pathway in file_pathways:
                pathway_scores[pathway].update(file_pathways[pathway])
        os.remove(file)

    return pathway_scores


def knockout_folder_pathway_propagations(knockout_folder: str, output_folder: str):
    files = sorted(list(Path(knockout_folder).glob("*.json")))
    folders = sorted(list(Path(knockout_folder).glob("*_randomized")))
    for i in range(len(files)):
        print(f"now handling {files[i]}")
        original_pathway_propagation = get_pathway_propagation_scores(files[i], folders[i],
                                                                      r"C:\studies\code\NetProp\data\c2.cp.v7.4.entrez.gmt")
        with open(os.path.join(output_folder, files[i].stem + "_pathways.json"), 'w') as handler:
            json.dump(original_pathway_propagation, handler, indent=4)


if __name__ == "__main__":
    knockout_folder_pathway_propagations(r'C:\studies\code\NetProp\src\temp\propagations\high_propagation_knockouts',
                                         r'C:\studies\code\NetProp\src\temp\knockout_pathways_enrichment\high_propagation')
    # pathway_propagations = \
    #     get_pathway_propagation_scores(r"C:\studies\code\NetProp\src\temp\propagations\knockouts\merged_covid_source.json",
    #                                    r"C:\studies\code\NetProp\src\temp\propagations\knockouts\merged_covid_source_randomized",
    #                                    r"C:\studies\code\NetProp\data\c2.cp.v7.4.entrez.gmt")
    # assert all([len(v) == 100 for v in pathway_propagations.values()])
    # with open(os.path.join(r"C:\studies\code\NetProp\src\temp\knockout_pathways_enrichment", "merged_covid_source_pathways.json"), 'w') as handler:
    #     json.dump(pathway_propagations, handler, indent=4)

