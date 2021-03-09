import json
from os import path

import ndex2.client
import networkx as nx
from pandas import read_csv
from mygene import MyGeneInfo

from datastructures.general_datastructures import HumanPriorMetadata, CovPriorMetadata
from datastructures.graph_datastructures import Protein

DEFAULT_COV_PROTEIN_ROLES_FILE_PATH = r"../data/cov_protein_roles.csv"
DEFAULT_HUMAN_PPI_FILE_PATH = r'../data/H_sapiens.net'


def acquire_cov_human_ppi(as_nx=True):
    ndex_server_url = 'public.ndexbio.org'
    hek293t_sars_cov_2_network_id =  "43803262-6d69-11ea-bfdc-0ac135e8bacf"

    nice_cx_from_server = ndex2.create_nice_cx_from_server(server=ndex_server_url,
                                                           uuid=hek293t_sars_cov_2_network_id)
    return nice_cx_from_server.to_networkx(mode="default") if as_nx else nice_cx_from_server


def init_cov_protein_roles(file_path):
    csv_dataframe = read_csv(file_path, index_col="protein")
    return {protein.lower(): set(column for column in csv_dataframe.columns if not roles.isnull()[column])
            for protein, roles in csv_dataframe.iterrows()}


# cov_human_ppi is a ppi where all edges originate in a cov protein and target a human protein
def prior_set_data_from_cov_human_ppi(cov_human_ppi, cov_protein_roles, from_targets=False):
    # In this PPI a gene is human if and only if it has an out degree of 0
    prior_set_data = dict()
    cov_human_interactions = [edge for edge in cov_human_ppi.edges(data=True)]
    for interaction in cov_human_interactions:
        interaction_name = interaction[2]["name"].lower().split()
        source, target = interaction_name[0], interaction_name[1]

        prior = target if from_targets else source
        metadata_class = HumanPriorMetadata if from_targets else CovPriorMetadata
        if prior not in prior_set_data:
            prior_set_data[prior] = metadata_class(name=prior)
        prior_set_data[prior].add_infection_roles(cov_protein_roles.get(source, set()))
        if from_targets:
            prior_set_data[prior].add_sources(source)
        else:
            prior_set_data[prior].add_targets(target)


    # map the i.d of each protein to its data. for humans proteins use entrezgene ids. for cov proteins arbitrary negative numbers (to differentiate from human proteins)
    gene_names = list(prior_set_data.keys())
    if from_targets:
        ncbi_query = MyGeneInfo().querymany(gene_names, scopes="symbol", fields=["entrezgene", "symbol"], species="human")
        for result in ncbi_query:
            prior_set_data[result["entrezgene"]] = prior_set_data.pop(result["symbol"].lower())
    else:
        for i in range(len(gene_names)):
            prior_set_data[-1*i - 1] = prior_set_data.pop(gene_names[i])

    return prior_set_data


def setup_cov_prior_set(human_ppi, cov_human_ppi, cov_protein_roles):
    """
    adds SARS-COV-19 proteins to a human ppi with MIST edges, and returns the set of these proteins as a prior set
    :param human_ppi: networkx graph where each node is a Protein with entrezgene id and edges are weighted by MIST
    :param cov_human_ppi: netowrkx graph taken from source cited in README.txt
    :param cov_protein_roles: dictionary that maps cov proteins to their roles in infecting human cells.
    :return:
    """
    interaction_dict = dict()
    prior_set = list()
    human_proteins_map = {protein.id: protein for protein in human_ppi}
    cov_human_interactions = [edge for edge in cov_human_ppi.edges(data=True)]
    for interaction in cov_human_interactions:
        interaction_name = interaction[2]["name"].lower().split()
        source, target = interaction_name[0], interaction_name[-1]
        mist_score = interaction[2].get("MIST", 0)
        if float(mist_score) > 0:
            if source not in interaction_dict:
                interaction_dict[source] = [(target, mist_score)]
            else:
                interaction_dict[source].append((target, mist_score))

    sources = list(interaction_dict.keys())
    with open(r'../data/gene_symbol_to_entrezgeneid.json', 'r') as json_handler:
        human_gene_symbol_to_id = json.load(json_handler)
    for i in range(len(sources)):
        source_name = sources[i]
        source_id = -1*i - 1
        source_roles = cov_protein_roles.get(source_name, set())
        # network_node = Protein(id=source_id, source_of=source_roles)
        network_node = Protein(id=source_id, source_of={source_id})
        prior_set.append(network_node)

        human_ppi.add_node(network_node)
        for interaction in interaction_dict[source_name]:
            target_symbol = interaction[0]
            target_id = int(human_gene_symbol_to_id[target_symbol])
            mist_score = float(interaction[1])
            human_ppi.add_edge(network_node, human_proteins_map[target_id], weight=mist_score)

    return prior_set


def acquire_prior_set_data(from_json=False, cov_protein_roles_path=r"../data/cov_protein_roles.csv", prior_set_from_human_targets=False):
    assert path.isfile(cov_protein_roles_path), "cov_protein_roles must be set to the path of an appropriate csv file"
    cov_human_ppi = acquire_cov_human_ppi()
    cov_protein_roles = init_cov_protein_roles(cov_protein_roles_path)
    if prior_set_from_human_targets:
        return prior_set_data_from_cov_human_ppi(cov_human_ppi, cov_protein_roles, from_targets=prior_set_from_human_targets)
    else:
        raise NotImplemented


def graph_from_file(file_path):
    source_node_index = 0
    target_node_index = 1
    edge_weight_index = 2

    graph = nx.DiGraph()
    discovered_nodes_map = dict()
    with open(file_path, 'r') as handler:
        for line in handler.readlines():
            values = line.split()
            for node_index in (values[source_node_index], values[target_node_index]):
                if node_index not in discovered_nodes_map:
                    discovered_nodes_map[node_index] = Protein(id=int(node_index))

            graph.add_edge(discovered_nodes_map[values[source_node_index]],
                           discovered_nodes_map[values[target_node_index]],
                           weight=float(values[edge_weight_index]))
    return graph



