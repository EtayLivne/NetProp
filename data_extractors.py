import ndex2.client
import networkx as nx
from mygene import MyGeneInfo
import pandas as pd
HEK293T_SARS_CoV_2__Network_id = "43803262-6d69-11ea-bfdc-0ac135e8bacf"
KROGAN_NODE_ID = 0
KROGAN_NODE_DATA = 1
NCBI_NODE_ID = 1

"""
anon_ndex=ndex2.client.Ndex2("http://public.ndexbio.org")

#cov_network = anon_ndex.get_network_as_cx_stream(HEK293T_SARS_CoV_2__Network_id)
nice_cx_from_server = ndex2.create_nice_cx_from_server(server='public.ndexbio.org', uuid=HEK293T_SARS_CoV_2__Network_id)

cov_networkx = nice_cx_from_server.to_networkx(mode="default")



"""

def acquire_cov_human_ppi(as_nx=True):
    nice_cx_from_server = ndex2.create_nice_cx_from_server(server='public.ndexbio.org',
                                                           uuid=HEK293T_SARS_CoV_2__Network_id)
    return nice_cx_from_server.to_networkx(mode="default") if as_nx else nice_cx_from_server

# cov_human_ppi is a ppi where all edges originate in a cov protein and target a human protein
def acquire_cov_human_ppi_prior_set_ids(cov_human_ppi):
    # In this PPI a gene is human if and only if it has an out degree of 0
    gene_names = list(node[KROGAN_NODE_DATA]["name"] for node in cov_human_ppi.nodes(data=True)
                      if cov_human_ppi.out_degree[node[KROGAN_NODE_ID]] == 0)
    ncbi_query = pd.DataFrame(MyGeneInfo().querymany(gene_names, scopes="symbol", fields="entrezgene", species="human"))
    return set(int(node_data[NCBI_NODE_ID]) for node_data in ncbi_query.values)


def acquire_cov_huamn_ppi_target_set(cov_human_ppi):