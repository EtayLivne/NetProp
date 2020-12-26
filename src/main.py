import json
from math import sqrt, pow

from numpy import sign

from data_extractors import *
from propagater import Propagater
from datastructures.graph_datastructures import is_cov_protein




# If node in 2 has less liquid than same node in 1, then the distance grew, and will be positive
def measure_propagation_distances(propagater1, propagater2, l2_distance=True):
    distance_dict = {}
    if l2_distance:
        def distance_function(x, y):
            return sign(x-y)*pow((x-y), 2)
    else:
        def distance_function(x, y):
            return x - y

    propagater2_nodes = {node.id: node for node in propagater2.graph.nodes}
    for node in propagater1.graph.nodes:
        if node.id in propagater2_nodes:
            for liquid_type, liquid in node.liquids.items():
                distance_dict[liquid_type] = distance_dict.get(liquid_type, 0) + \
                                             distance_function(liquid, propagater2_nodes[node.id].liquids.get(liquid_type, 0))

    for liquid_type in distance_dict:
        distance_dict[liquid_type] = sign(distance_dict[liquid_type]) * sqrt(abs(distance_dict[liquid_type]))

    return distance_dict


def measure_gene_knockout_impact(graph, knockout_gene_ids, confidence_coef,
                                 halt_condition_gap=Propagater._DEFUALT_HALT_CONDITION_GAP, prior_set_ids=None):
    if prior_set_ids is None:
        prior_set = {node for node in graph.nodes if is_cov_protein(node)}
    else:
        protein_set = {Protein(id=prior_id) for prior_id in prior_set_ids}

    propagater = Propagater(graph, confidence_coef, halt_condition_gap=halt_condition_gap)
    print("propagating original graph")
    propagation_distance_dict_l1 = dict()
    propagation_distance_dict_l2 = dict()
    propagater.propagate(prior_set)
    print("done propagating original graph")
    print(f'There are {len(knockout_gene_ids)} knockouts to measure')
    i = 0
    for knockout_genes in knockout_gene_ids:
        i = i + 1
        knockout_propagater = propagater.gene_knockout_copy(knockout_genes)
        print(f'begininig propagation on network that knocked out {knockout_genes}')
        knockout_propagater.propagate(prior_set)
        print(f'done with knockout {i} out of {len(knockout_gene_ids)}')
        propagation_distance_dict_l1[str(knockout_genes)] = \
            measure_propagation_distances(propagater, knockout_propagater, l2_distance=False)
        propagation_distance_dict_l2[str(knockout_genes)] = \
            measure_propagation_distances(propagater, knockout_propagater)

    return propagation_distance_dict_l1, propagation_distance_dict_l2


def prioritize_knockout_gene_sets(propagation_distance_dict, l2_distance=True):
    if l2_distance:
        def norm(distances):
            sum_with_signs = sum(sign(d) * pow(d, 2) for d in distances)

            return sign(sum_with_signs) * sqrt(abs(sum_with_signs))
    else:
        def norm(distances):
            return sum(d for d in distances)

    return sorted([(knockout_set, norm(distances.values()), distances) for knockout_set, distances in propagation_distance_dict.items()],
                  key=lambda item: item[1], reverse=True)


def compare_priorities(*propagation_distance_dicts, top_k=100):
    rankings_dict = dict()
    for d in propagation_distance_dicts:
        distance_dict_name = d[0]
        distance_dict = d[1]
        l1_ranking = prioritize_knockout_gene_sets(distance_dict, l2_distance=False)
        l2_ranking = prioritize_knockout_gene_sets(distance_dict)
        rankings_dict[distance_dict_name] = {
            f'top_{top_k}_l1_distances': [{
                    f'norm': e[1],
                    f'distances': e[2],
                    "knocked_genes_entrezids": e[0],
                    "knockout_gene_ids": [int(num) for num in e[0][1:-1].split(',')]
            } for e in l1_ranking[:top_k]],
            f'top_{top_k}_l2_distances': [{
                    f'norm': e[1],
                    f'distances': e[2],
                    "knocked_genes_entrezids": e[0],
                    "knockout_gene_ids": [int(num) for num in e[0][1:-1].split(',')]
            } for e in l2_ranking[:top_k]],
        }

    return rankings_dict



if __name__ == "__main__":
    import json
    with open(r'D:\Etay\studies\thesis\code\NetProp\data\COVID-19-–-VeroE6-IF_reframedb-data_2020-12-13.json', 'r') as h:
        d = json.load(h)
    with open(r'D:\Etay\studies\thesis\code\NetProp\data\COVID-19-–-VeroE6-IF_reframedb-data_2020-12-13.json', 'w') as h:
        json.dump(d, h, indent=4)
    # cov_human_ppi = acquire_cov_human_ppi()
    # cov_protein_roles = init_cov_protein_roles(DEFAULT_COV_PROTEIN_ROLES_FILE_PATH)
    # human_ppi = graph_from_file(r"..\data\H_sapiens.net")
    # print("graph created!")
    # prior_set = setup_cov_prior_set(human_ppi, cov_human_ppi, cov_protein_roles)
    #
    #
    # prop = Propagater(human_ppi, 0.9)
    # prop.propagate(prior_set)
    # p_list = sorted(list(node for node in prop.graph.nodes if not is_cov_protein(node)),
    #                 key=lambda node: sum(node.liquids.values()), reverse=True)[:50]
    # # knockout_gene_ids = [{p.id} for p in p_list]
    # knockout_gene_ids = []
    # for i in range(len(p_list)):
    #     knockout_gene_ids.extend([{p_list[i].id, p_list[j].id} for j in range(i + 1, len(p_list))])
    # gene_knockout_distance_dicts = measure_gene_knockout_impact(human_ppi, knockout_gene_ids, 0.9)
    # rankings = compare_priorities(["l1_distance_dict", gene_knockout_distance_dicts[0]],
    #                               ["l2_distance_dict", gene_knockout_distance_dicts[1]])
    #
    # with open(r'../results/l1_l2_double_knockout_signed_l2_norm.json', "w") as handler:
    #     json.dump(rankings, handler, indent=4)
    #
    # jhgf  = 9