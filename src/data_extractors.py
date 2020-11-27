import ndex2.client

from os import path
from pandas import read_csv
from mygene import MyGeneInfo

from datastructures.general_datastructures import CovTargetMetadata



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
def acquire_cov_human_ppi_prior_set_data(cov_human_ppi, cov_protein_roles):
    # In this PPI a gene is human if and only if it has an out degree of 0
    prior_set_data = dict()
    cov_human_interactions = [edge for edge in cov_human_ppi.edges(data=True)]
    for interaction in cov_human_interactions:
        interaction_name = interaction[2]["name"].lower().split()
        source = interaction_name[0]
        target = interaction_name[-1]
        if target not in prior_set_data:
            prior_set_data[target] = CovTargetMetadata(name=target)

        prior_set_data[target].add_sources(source)
        prior_set_data[target].add_infection_roles(cov_protein_roles.get(source, set()))

    human_gene_names = list(prior_set_data.keys())
    ncbi_query = MyGeneInfo().querymany(human_gene_names, scopes="symbol", fields=["entrezgene", "symbol"], species="human")
    for result in ncbi_query:
        prior_set_data[result["entrezgene"]] = prior_set_data.pop(result["symbol"].lower())

    return prior_set_data


def acquire_prior_set_data(from_json=False, cov_protein_roles_path=r"../data/cov_protein_roles.csv"):
    assert(path.isfile(cov_protein_roles_path), "cov_protein_roles must be set to the path of an appropriate csv file")
    cov_human_ppi = acquire_cov_human_ppi()
    cov_protein_roles = init_cov_protein_roles(cov_protein_roles_path)
    return acquire_cov_human_ppi_prior_set_data(cov_human_ppi, cov_protein_roles)


