import json
from math import sqrt, pow
import re
from multiprocessing import Process
from numpy import sign
from random import choice
import json
import os
from data_extractors import *
from propagater import Propagater, PropagationContainer, PropagationGraph
from datastructures.graph_datastructures import is_cov_protein
from datastructures.drugbank_data import DrugBankProteinTargetsData
# from datastructures.general_datastructures import KnockoutGeneSet, SymbolEntrezgeneMap, init_symbol_entrezgene_map
from common.data_handlers.managers import HSapiensManager
from common.data_handlers.extractors import GeneInfoExtractor





# order of propagaters matter - increases in liquid from first to second will be counted as positive, decreases as negative
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


def measure_gene_knockout_impact(graph, knockout_gene_sets, confidence_coef,
                                 halt_condition_gap=Propagater._DEFAULT_HALT_CONDITION_GAP, prior_set_ids=None):
    if prior_set_ids is None:
        prior_set = {node for node in graph.nodes if is_cov_protein(node)}
    else:
        prior_set = {Protein(id=prior_id) for prior_id in prior_set_ids}

    propagater = Propagater(graph, confidence_coef, halt_condition_gap=halt_condition_gap)
    print("propagating original graph")
    propagation_distance_dict_l1 = dict()
    propagation_distance_dict_l2 = dict()
    propagater.propagate(prior_set)
    print("done propagating original graph")
    print(f'There are {len(knockout_gene_sets)} knockouts to measure')
    i = 0
    for knockout_gene_set in knockout_gene_sets:
        i = i + 1
        knockout_propagater = propagater.gene_knockout_copy(knockout_gene_set.gene_set)
        print(f'beginning propagation on network that knocked out {knockout_gene_set.name}')
        knockout_propagater.propagate(prior_set)
        print(f'done with knockout {i} out of {len(knockout_gene_sets)}')
        propagation_distance_dict_l1[str(knockout_gene_set.name)] = \
            measure_propagation_distances(propagater, knockout_propagater, l2_distance=False)
        propagation_distance_dict_l2[str(knockout_gene_set.name)] = \
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


def compare_priorities(knockout_gene_sets, *propagation_distance_dicts, top_k=100):

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
                    "knockout_gene_ids": list(next(gs for gs in knockout_gene_sets if gs.name==e[0]).gene_set)
            } for e in l1_ranking[:top_k]],
            f'top_{top_k}_l2_distances': [{
                    f'norm': e[1],
                    f'distances': e[2],
                    "knocked_genes_entrezids": e[0],
                    "knockout_gene_ids": list(next(gs for gs in knockout_gene_sets if gs.name==e[0]).gene_set)
            } for e in l2_ranking[:top_k]],
        }

    return rankings_dict


def generate_knockout_profiles_from_drugs():
    drugbank_data = DrugBankProteinTargetsData()
    drugbank_data.init_from_file(r'../data/all_drugbank_drug_targets.csv')
    drug_to_target_map = drugbank_data.get_human_drug_protein_targets()
    symbol_entrezgene_map = SymbolEntrezgeneMap()
    profiles_list = []
    for drug, targets in drug_to_target_map.items():
        gene_set = set()
        for target in targets:
            try:
                gene_set.add(symbol_entrezgene_map.get_entrezgene(target))
            except KeyError:
                pass
        profiles_list.append(KnockoutGeneSet(name=drug, gene_set=gene_set))

    return profiles_list


def generate_rank_equivalent_random_network(network):
    edges_to_switch = 10 * len(list(network.edges))
    switch_types = ["place", "target", "source"]
    edge_source_sample = []
    for node in network.nodes:
        edge_source_sample.extend([node] * network.degree(node))
    while edges_to_switch:
        source_1 = choice(edge_source_sample)
        target_1 = choice(list(network.neighbors(source_1)))
        random_edge_1 = (source_1, target_1)
        random_edge_1_attrs = dict(network.edges[random_edge_1])
        while True:
            source_2 = choice(edge_source_sample)
            target_2 = choice(list(network.neighbors(source_2)))
            if source_2 not in random_edge_1 and target_2 not in  random_edge_1:
                break
        random_edge_2 = (source_2, target_2)
        random_edge_2_attrs = dict(network.edges[random_edge_2])

        switch_type = choice(switch_types)
        if switch_type == "place":
            network.remove_edge(*random_edge_1)
            network.remove_edge(*random_edge_2)
            network.add_edge(*random_edge_1, **random_edge_2_attrs)
            network.add_edge(*random_edge_2, **random_edge_1_attrs)
        if switch_type == "target":
            if network.has_edge(random_edge_1[0], random_edge_2[1]) or network.has_edge(random_edge_2[0], random_edge_1[1]):
                continue
            network.remove_edge(*random_edge_1)
            network.remove_edge(*random_edge_2)
            network.add_edge(random_edge_1[0], random_edge_2[1], **random_edge_1_attrs)
            network.add_edge(random_edge_2[0], random_edge_1[1], **random_edge_2_attrs)
        if switch_type == "source":
            if network.has_edge(random_edge_2[0], random_edge_1[1]) or network.has_edge(random_edge_1[0], random_edge_2[1]):
                continue
            network.remove_edge(*random_edge_1)
            network.remove_edge(*random_edge_2)
            network.add_edge(random_edge_2[0], random_edge_1[1], **random_edge_1_attrs)
            network.add_edge(random_edge_1[0], random_edge_2[1], **random_edge_2_attrs)

        edges_to_switch -= 1
        if edges_to_switch % 1000 == 0:
            print(f'{edges_to_switch} edges remain to switch')


def randomize_h_sapiens():
    h_sapiens = HSapiensManager(r"C:\studies\code\NetProp\data\H_sapiens.net")
    human_ppi = h_sapiens.get_data()
    generate_rank_equivalent_random_network(human_ppi)
    return human_ppi


def prior_set_from_json(path_to_json, prior_set_identifier):
    with open(path_to_json, 'r') as handler:
        data = json.load(handler)

    try:
        gene_names = data[prior_set_identifier]
    except KeyError:
        raise KeyError(f"no prior set with identifier {prior_set_identifier} in json file {path_to_json}")

    gene_ids_extractor = GeneInfoExtractor()
    query_results = [result for result in gene_ids_extractor.extract(gene_names) if "entrezgene" in result]
    return [int(result['entrezgene']) for result in query_results]


def propagate_on_random_networks(variant_name, with_ph):
    ph_modified_name = "ph_" + variant_name if with_ph else variant_name
    variant_prior_sets_file_name = "ph_variant_prior_sets" if with_ph else "variant_prior_sets"
    output_dir = r'C:\studies\code\NetProp\results\variants_blitz\randomized_networks_propagations\{}'.format(ph_modified_name)
    data_dir = r'C:\studies\code\NetProp\results\variants_blitz\randomized_networks'
    variant_prior_set = prior_set_from_json(r"C:\studies\code\NetProp\data\variants_blitz\{}.json".format(variant_prior_sets_file_name),
                                            f"{variant_name}")
    for randomized_human_ppi_file in os.listdir(data_dir):
        num = (re.findall("\d+", randomized_human_ppi_file))[-1]
        if not num:
            continue
        if int(num) < 74:
            continue
        h_sapiens_file_path = data_dir + f"\\{randomized_human_ppi_file}"
        h_sapiens = HSapiensManager(h_sapiens_file_path)
        human_ppi = h_sapiens.get_data()
        problems = []
        for gene_id in variant_prior_set:
            try:
                human_ppi.nodes[gene_id][PropagationGraph.CONTAINER_PROPERTIES_KEY].source_of = {"liquid"}
            except KeyError:
                print(f"gene {gene_id} for prior set does not, in fact, appear to be a human gene")
                problems.append(gene_id)
        print(f"ignored the following genes: {problems}")
        variant_prior_set = [gene_id for gene_id in variant_prior_set if gene_id not in problems]

        prop = Propagater(human_ppi, 0.1)
        prop.propagate(variant_prior_set)
        node_tuples = sorted([(node, data[prop._LIQUIDS]["liquid"]) for node, data in prop.graph.nodes(data=True)],
                             key=lambda x: x[1], reverse=True)
        output_path = output_dir + f"\\{ph_modified_name}_{randomized_human_ppi_file}.json"
        with open(output_path, 'w') as handler:
            json.dump([{"gene_id": node_tuple[0], "liquid": node_tuple[1], "in_prior_set": node_tuple[0] in variant_prior_set} for node_tuple in node_tuples], handler, indent=4)
        print("yo")


def generic_propagate_on_random_networks(prior_set, data_dir, output_dir, output_file_name_prefix):

    for randomized_human_ppi_file in os.listdir(data_dir):
        h_sapiens_file_path = data_dir + f"\\{randomized_human_ppi_file}"
        h_sapiens = HSapiensManager(h_sapiens_file_path)
        human_ppi = h_sapiens.get_data()
        problems = []
        for gene_id in prior_set:
            try:
                human_ppi.nodes[gene_id][PropagationGraph.CONTAINER_PROPERTIES_KEY].source_of = {"liquid"}
            except KeyError:
                print(f"gene {gene_id} for prior set does not, in fact, appear to be a human gene")
                problems.append(gene_id)
        print(f"ignored the following genes: {problems}")
        prior_set = [gene_id for gene_id in prior_set if gene_id not in problems]

        prop = Propagater(human_ppi, 0.1)
        prop.propagate(prior_set)
        node_tuples = sorted([(node, data[prop._LIQUIDS]["liquid"]) for node, data in prop.graph.nodes(data=True)],
                             key=lambda x: x[1], reverse=True)
        output_path = output_dir + f"\\{output_file_name_prefix}_{randomized_human_ppi_file}.json"
        with open(output_path, 'w') as handler:
            json.dump([{"gene_id": node_tuple[0], "liquid": node_tuple[1], "in_prior_set": node_tuple[0] in prior_set} for node_tuple in node_tuples], handler, indent=4)
        print("yo")


if __name__ == "__main__":
    # init_symbol_entrezgene_map(from_file=r'../data/gene_symbol_to_entrezgeneid_new.json')
    # # with open(r'D:\Etay\studies\thesis\code\NetProp\data\COVID-19-–-VeroE6-IF_reframedb-data_2020-12-13.json', 'r') as h:
    # #     d = json.load(h)
    # # with open(r'D:\Etay\studies\thesis\code\NetProp\data\COVID-19-–-VeroE6-IF_reframedb-data_2020-12-13.json', 'w') as h:
    # #     json.dump(d, h, indent=4)
    # cov_human_ppi = acquire_cov_human_ppi()
    # cov_protein_roles = init_cov_protein_roles(DEFAULT_COV_PROTEIN_ROLES_FILE_PATH)
    # human_ppi = graph_from_file(r"..\data\H_sapiens.net")
    # print("graph created!")
    # prior_set = setup_cov_prior_set(human_ppi, cov_human_ppi, cov_protein_roles)
    #
    #
    # prop = Propagater(human_ppi, 0.9)
    # prop.propagate(prior_set)
    # # p_list = sorted(list(node for node in prop.graph.nodes if not is_cov_protein(node)),
    # #                 key=lambda node: sum(node.liquids.values()), reverse=True)[:50]
    # # knockout_gene_ids = [{p.id} for p in p_list]
    # # knockout_gene_ids = []
    # # for i in range(len(p_list)):
    # #     knockout_gene_ids.extend([{p_list[i].id, p_list[j].id} for j in range(i + 1, len(p_list))])
    # knockout_gene_sets = generate_knockout_profiles_from_drugs()
    # for i in range(int(len(knockout_gene_sets)/100)):
    #     gene_knockout_distance_dicts = measure_gene_knockout_impact(human_ppi, knockout_gene_sets[100*i:100*(i+1)], 0.9)
    #     rankings = compare_priorities(knockout_gene_sets,
    #                                   ["l1_distance_dict", gene_knockout_distance_dicts[0]],
    #                                   ["l2_distance_dict", gene_knockout_distance_dicts[1]])
    #
    #     with open(r'../results/drugs {} to {}.json'.format(100*i, 100*(i+1)), "w") as handler:
    #         json.dump(rankings, handler, indent=4)


    ######### refactor monkey testing #########
    # h_sapiens = HSapiensManager(r"C:\studies\code\NetProp\data\H_sapiens.net")
    # human_ppi = h_sapiens.get_data()
    # human_ppi.nodes[1][PropagationGraph.CONTAINER_PROPERTIES_KEY].source_of = {"liquid"}
    # prior_set = [1]
    # prop = Propagater(human_ppi, 0.1)
    # prop.propagate(prior_set)
    # node_tuples = sorted([(node, data[prop._LIQUIDS]["liquid"]) for node, data in prop.graph.nodes(data=True)],
    #                      key=lambda x: x[1])
    # with open('output.txt', 'w') as handler:
    #         handler.write("\n".join(str(node_tuple) for node_tuple in node_tuples))
    # print("yo")

    ######## generating randomized networks for p value test #########
        # for i in range(26):
        #     randomized_network = randomize_h_sapiens()
        #     file_path = r'C:\studies\code\NetProp\results\variants_blitz\randomized_networks\randomized_h_sapiens_{}.net'.format(74 + i)
        #     with open(file_path, 'w') as handler:
        #         for edge in randomized_network.edges(data=True):
        #             handler.write(f'{edge[0]} {edge[1]} {edge[2]["weight"]} 0\n')
        #     print(f"\n\n********done writing file {i}********\n\n")
    




    ######### propagating on human ppi with actual variant differentiated genes as prior set #########
    # variant_prior_set = prior_set_from_json(r"C:\studies\code\NetProp\data\variants_blitz\ph_variant_prior_sets.json", "VIC")
    # h_sapiens = HSapiensManager(r"C:\studies\code\NetProp\data\H_sapiens.net")
    # human_ppi = h_sapiens.get_data()
    # problems = []
    # for gene_id in variant_prior_set:
    #     try:
    #         human_ppi.nodes[gene_id][PropagationGraph.CONTAINER_PROPERTIES_KEY].source_of = {"liquid"}
    #     except KeyError:
    #         print(f"gene {gene_id} for prior set does not, in fact, appear to be a human gene")
    #         problems.append(gene_id)
    # print(f"ignored the following genes: {problems}")
    # variant_prior_set = [gene_id for gene_id in variant_prior_set if gene_id not in problems]
    #
    # prop = Propagater(human_ppi, 0.1)
    # prop.propagate(variant_prior_set)
    # node_tuples = sorted([(node, data[prop._LIQUIDS]["liquid"]) for node, data in prop.graph.nodes(data=True)],
    #                      key=lambda x: x[1], reverse=True)
    # with open('ph_VIC_output.json', 'w') as handler:
    #     json.dump([{"gene_id": node_tuple[0], "liquid": node_tuple[1], "in_prior_set": node_tuple[0] in variant_prior_set} for node_tuple in node_tuples], handler, indent=4)
    # print("yo")


    ######### propagating on human ppi with combined variant differentiated genes as prior set #########
    # VIC_prior_set = prior_set_from_json(r"C:\studies\code\NetProp\data\variants_blitz\variant_prior_sets.json", "VIC")
    # Kent_prior_set = prior_set_from_json(r"C:\studies\code\NetProp\data\variants_blitz\variant_prior_sets.json",
    #                                     "Kent")
    # EU_prior_set = prior_set_from_json(r"C:\studies\code\NetProp\data\variants_blitz\variant_prior_sets.json",
    #                                     "EU120")
    # # Kent_minus_VIC = list(set(Kent_prior_set) - set(VIC_prior_set))
    # # Kent_and_EU_minus_VIC = list(set(Kent_prior_set) & set(EU_prior_set) - set(VIC_prior_set))
    # VIC_minus_Kent = list(set(VIC_prior_set) - set(Kent_prior_set))
    # VIC_minus_Kent_and_EU = list(set(VIC_prior_set) - set(Kent_prior_set) & set(EU_prior_set))
    # h_sapiens = HSapiensManager(r"C:\studies\code\NetProp\data\H_sapiens.net")
    # human_ppi = h_sapiens.get_data()
    # problems = []
    # for gene_id in VIC_minus_Kent:
    #     try:
    #         human_ppi.nodes[gene_id][PropagationGraph.CONTAINER_PROPERTIES_KEY].source_of = {"liquid"}
    #     except KeyError:
    #         print(f"gene {gene_id} for prior set does not, in fact, appear to be a human gene")
    #         problems.append(gene_id)
    # print(f"ignored the following genes: {problems}")
    # VIC_minus_Kent = [gene_id for gene_id in VIC_minus_Kent if gene_id not in problems]
    #
    # prop = Propagater(human_ppi, 0.1)
    # prop.propagate(VIC_minus_Kent)
    # node_tuples = sorted([(node, data[prop._LIQUIDS]["liquid"]) for node, data in prop.graph.nodes(data=True)],
    #                      key=lambda x: x[1], reverse=True)
    # with open('VIC_minus_Kent.json', 'w') as handler:
    #     json.dump([{"gene_id": node_tuple[0], "liquid": node_tuple[1], "in_prior_set": node_tuple[0] in VIC_minus_Kent} for node_tuple in node_tuples], handler, indent=4)
    # print("yo")
    # h_sapiens = HSapiensManager(r"C:\studies\code\NetProp\data\H_sapiens.net")
    # human_ppi = h_sapiens.get_data()
    # problems = []
    # for gene_id in VIC_minus_Kent_and_EU:
    #     try:
    #         human_ppi.nodes[gene_id][PropagationGraph.CONTAINER_PROPERTIES_KEY].source_of = {"liquid"}
    #     except KeyError:
    #         print(f"gene {gene_id} for prior set does not, in fact, appear to be a human gene")
    #         problems.append(gene_id)
    # print(f"ignored the following genes: {problems}")
    # VIC_minus_Kent_and_EU = [gene_id for gene_id in VIC_minus_Kent_and_EU if gene_id not in problems]
    #
    # prop = Propagater(human_ppi, 0.1)
    # prop.propagate(VIC_minus_Kent_and_EU)
    # node_tuples = sorted([(node, data[prop._LIQUIDS]["liquid"]) for node, data in prop.graph.nodes(data=True)],
    #                      key=lambda x: x[1], reverse=True)
    # with open('VIC_minus_Kent_and_EU.json', 'w') as handler:
    #     json.dump([{"gene_id": node_tuple[0], "liquid": node_tuple[1], "in_prior_set": node_tuple[0] in VIC_minus_Kent_and_EU} for node_tuple in node_tuples], handler, indent=4)
    # print("yo")
    #
    #
    #
    # ph_VIC_prior_set = prior_set_from_json(r"C:\studies\code\NetProp\data\variants_blitz\ph_variant_prior_sets.json", "VIC")
    # ph_Kent_prior_set = prior_set_from_json(r"C:\studies\code\NetProp\data\variants_blitz\ph_variant_prior_sets.json",
    #                                     "Kent")
    # ph_EU_prior_set = prior_set_from_json(r"C:\studies\code\NetProp\data\variants_blitz\ph_variant_prior_sets.json",
    #                                     "EU120")
    # ph_VIC_minus_Kent = list(set(ph_VIC_prior_set) - set(ph_Kent_prior_set))
    # ph_VIC_minus_Kent_and_EU = list(set(ph_VIC_prior_set) - set(ph_Kent_prior_set) & set(ph_EU_prior_set))
    # h_sapiens = HSapiensManager(r"C:\studies\code\NetProp\data\H_sapiens.net")
    # human_ppi = h_sapiens.get_data()
    # problems = []
    # for gene_id in ph_VIC_minus_Kent:
    #     try:
    #         human_ppi.nodes[gene_id][PropagationGraph.CONTAINER_PROPERTIES_KEY].source_of = {"liquid"}
    #     except KeyError:
    #         print(f"gene {gene_id} for prior set does not, in fact, appear to be a human gene")
    #         problems.append(gene_id)
    # print(f"ignored the following genes: {problems}")
    # ph_VIC_minus_Kent = [gene_id for gene_id in ph_VIC_minus_Kent if gene_id not in problems]
    #
    # prop = Propagater(human_ppi, 0.1)
    # prop.propagate(ph_VIC_minus_Kent)
    # node_tuples = sorted([(node, data[prop._LIQUIDS]["liquid"]) for node, data in prop.graph.nodes(data=True)],
    #                      key=lambda x: x[1], reverse=True)
    # with open('ph_VIC_minus_Kent.json', 'w') as handler:
    #     json.dump([{"gene_id": node_tuple[0], "liquid": node_tuple[1], "in_prior_set": node_tuple[0] in ph_VIC_minus_Kent} for node_tuple in node_tuples], handler, indent=4)
    # print("yo")
    # h_sapiens = HSapiensManager(r"C:\studies\code\NetProp\data\H_sapiens.net")
    # human_ppi = h_sapiens.get_data()
    # problems = []
    # for gene_id in ph_VIC_minus_Kent_and_EU:
    #     try:
    #         human_ppi.nodes[gene_id][PropagationGraph.CONTAINER_PROPERTIES_KEY].source_of = {"liquid"}
    #     except KeyError:
    #         print(f"gene {gene_id} for prior set does not, in fact, appear to be a human gene")
    #         problems.append(gene_id)
    # print(f"ignored the following genes: {problems}")
    # ph_VIC_minus_Kent_and_EU = [gene_id for gene_id in ph_VIC_minus_Kent_and_EU if gene_id not in problems]
    #
    # prop = Propagater(human_ppi, 0.1)
    # prop.propagate(ph_VIC_minus_Kent_and_EU)
    # node_tuples = sorted([(node, data[prop._LIQUIDS]["liquid"]) for node, data in prop.graph.nodes(data=True)],
    #                      key=lambda x: x[1], reverse=True)
    # with open('ph_VIC_minus_Kent_and_EU.json', 'w') as handler:
    #     json.dump([{"gene_id": node_tuple[0], "liquid": node_tuple[1], "in_prior_set": node_tuple[0] in ph_VIC_minus_Kent_and_EU} for node_tuple in node_tuples], handler, indent=4)
    # print("yo")




    ######### propagating on human ppi with actual variant differentiated genes as prior set, on randomized networks #########
    # output_dir = r'C:\studies\code\NetProp\results\variants_blitz\randomized_networks_propagations\ph_EU120'
    # data_dir = r'C:\studies\code\NetProp\results\variants_blitz\randomized_networks'
    # variant_prior_set = prior_set_from_json(r"C:\studies\code\NetProp\data\variants_blitz\variant_prior_sets.json",
    #                                         "EU120A")
    # for randomized_human_ppi_file in os.listdir(data_dir):
    #     if "super" not in randomized_human_ppi_file:
    #         continue
    #     h_sapiens_file_path = data_dir + f"\\{randomized_human_ppi_file}"
    #     h_sapiens = HSapiensManager(h_sapiens_file_path)
    #     human_ppi = h_sapiens.get_data()
    #     problems = []
    #     for gene_id in variant_prior_set:
    #         try:
    #             human_ppi.nodes[gene_id][PropagationGraph.CONTAINER_PROPERTIES_KEY].source_of = {"liquid"}
    #         except KeyError:
    #             print(f"gene {gene_id} for prior set does not, in fact, appear to be a human gene")
    #             problems.append(gene_id)
    #     print(f"ignored the following genes: {problems}")
    #     variant_prior_set = [gene_id for gene_id in variant_prior_set if gene_id not in problems]
    #
    #     prop = Propagater(human_ppi, 0.1)
    #     prop.propagate(variant_prior_set)
    #     node_tuples = sorted([(node, data[prop._LIQUIDS]["liquid"]) for node, data in prop.graph.nodes(data=True)],
    #                          key=lambda x: x[1], reverse=True)
    #     output_path = output_dir + f"\\ph_EU120_{randomized_human_ppi_file}.json"
    #     with open(output_path, 'w') as handler:
    #         json.dump([{"gene_id": node_tuple[0], "liquid": node_tuple[1], "in_prior_set": node_tuple[0] in variant_prior_set} for node_tuple in node_tuples], handler, indent=4)
    #     print("yo")

    # output_dir = r'C:\studies\code\NetProp\results\variants_blitz\randomized_networks_propagations\ph_Kent'
    # data_dir = r'C:\studies\code\NetProp\results\variants_blitz\randomized_networks'
    # variant_prior_set = prior_set_from_json(r"C:\studies\code\NetProp\data\variants_blitz\ph_variant_prior_sets.json",
    #                                         "Kent")
    # for randomized_human_ppi_file in os.listdir(data_dir):
    #     num = (re.findall("\d+", randomized_human_ppi_file))[0]
    #     if not num or int(num) < 74:
    #         continue
    #     h_sapiens_file_path = data_dir + f"\\{randomized_human_ppi_file}"
    #     h_sapiens = HSapiensManager(h_sapiens_file_path)
    #     human_ppi = h_sapiens.get_data()
    #     problems = []
    #     for gene_id in variant_prior_set:
    #         try:
    #             human_ppi.nodes[gene_id][PropagationGraph.CONTAINER_PROPERTIES_KEY].source_of = {"liquid"}
    #         except KeyError:
    #             print(f"gene {gene_id} for prior set does not, in fact, appear to be a human gene")
    #             problems.append(gene_id)
    #     print(f"ignored the following genes: {problems}")
    #     variant_prior_set = [gene_id for gene_id in variant_prior_set if gene_id not in problems]
    #
    #     prop = Propagater(human_ppi, 0.1)
    #     prop.propagate(variant_prior_set)
    #     node_tuples = sorted([(node, data[prop._LIQUIDS]["liquid"]) for node, data in prop.graph.nodes(data=True)],
    #                          key=lambda x: x[1], reverse=True)
    #     output_path = output_dir + f"\\ph_Kent_{randomized_human_ppi_file}.json"
    #     with open(output_path, 'w') as handler:
    #         json.dump([{"gene_id": node_tuple[0], "liquid": node_tuple[1],
    #                     "in_prior_set": node_tuple[0] in variant_prior_set} for node_tuple in node_tuples], handler,
    #                   indent=4)
    #     print("yo")
    #
    # output_dir = r'C:\studies\code\NetProp\results\variants_blitz\randomized_networks_propagations\ph_VIC'
    # data_dir = r'C:\studies\code\NetProp\results\variants_blitz\randomized_networks'
    # variant_prior_set = prior_set_from_json(
    #     r"C:\studies\code\NetProp\data\variants_blitz\ph_variant_prior_sets.json",
    #     "VIC")
    # for randomized_human_ppi_file in os.listdir(data_dir):
    #     num = (re.findall("\d+", randomized_human_ppi_file))[0]
    #     if not num or int(num) < 74:
    #         continue
    #     h_sapiens_file_path = data_dir + f"\\{randomized_human_ppi_file}"
    #     h_sapiens = HSapiensManager(h_sapiens_file_path)
    #     human_ppi = h_sapiens.get_data()
    #     problems = []
    #     for gene_id in variant_prior_set:
    #         try:
    #             human_ppi.nodes[gene_id][PropagationGraph.CONTAINER_PROPERTIES_KEY].source_of = {"liquid"}
    #         except KeyError:
    #             print(f"gene {gene_id} for prior set does not, in fact, appear to be a human gene")
    #             problems.append(gene_id)
    #     print(f"ignored the following genes: {problems}")
    #     variant_prior_set = [gene_id for gene_id in variant_prior_set if gene_id not in problems]
    #
    #     prop = Propagater(human_ppi, 0.1)
    #     prop.propagate(variant_prior_set)
    #     node_tuples = sorted(
    #         [(node, data[prop._LIQUIDS]["liquid"]) for node, data in prop.graph.nodes(data=True)],
    #         key=lambda x: x[1], reverse=True)
    #     output_path = output_dir + f"\\ph_VIC_{randomized_human_ppi_file}.json"
    #     with open(output_path, 'w') as handler:
    #         json.dump([{"gene_id": node_tuple[0], "liquid": node_tuple[1],
    #                     "in_prior_set": node_tuple[0] in variant_prior_set} for node_tuple in node_tuples],
    #                   handler, indent=4)
    #     print("yo")

    ######## multiprocessed propagating on human ppi with actual variant differentiated genes as prior set, on randomized networks #########

    VIC_prior_set = prior_set_from_json(r"C:\studies\code\NetProp\data\variants_blitz\variant_prior_sets.json", "VIC")
    Kent_prior_set = prior_set_from_json(r"C:\studies\code\NetProp\data\variants_blitz\variant_prior_sets.json",
                                        "Kent")
    EU_prior_set = prior_set_from_json(r"C:\studies\code\NetProp\data\variants_blitz\variant_prior_sets.json",
                                        "EU120")
    VIC_minus_Kent = list(set(VIC_prior_set) - set(Kent_prior_set))
    VIC_minus_Kent_and_EU = list(set(VIC_prior_set) - set(Kent_prior_set) & set(EU_prior_set))
    problems = []
    h_sapiens = HSapiensManager(r"C:\studies\code\NetProp\data\H_sapiens.net")
    human_ppi = h_sapiens.get_data()
    for gene_id in VIC_minus_Kent:
        try:
            human_ppi.nodes[gene_id][PropagationGraph.CONTAINER_PROPERTIES_KEY].source_of = {"liquid"}
        except KeyError:
            print(f"gene {gene_id} for prior set does not, in fact, appear to be a human gene")
            problems.append(gene_id)
    print(f"ignored the following genes: {problems}")
    VIC_minus_Kent = [gene_id for gene_id in VIC_minus_Kent if gene_id not in problems]
    for gene_id in VIC_minus_Kent_and_EU:
        try:
            human_ppi.nodes[gene_id][PropagationGraph.CONTAINER_PROPERTIES_KEY].source_of = {"liquid"}
        except KeyError:
            print(f"gene {gene_id} for prior set does not, in fact, appear to be a human gene")
            problems.append(gene_id)
    print(f"ignored the following genes: {problems}")
    VIC_minus_Kent_and_EU = [gene_id for gene_id in VIC_minus_Kent_and_EU if gene_id not in problems]


    ph_Kent_prior_set = prior_set_from_json(r"C:\studies\code\NetProp\data\variants_blitz\ph_variant_prior_sets.json",
                                        "Kent")
    ph_EU_prior_set = prior_set_from_json(r"C:\studies\code\NetProp\data\variants_blitz\ph_variant_prior_sets.json",
                                        "EU120")
    ph_VIC_prior_set = prior_set_from_json(r"C:\studies\code\NetProp\data\variants_blitz\ph_variant_prior_sets.json",
                                           "VIC")
    ph_VIC_minus_Kent = list(set(ph_Kent_prior_set) - set(ph_VIC_prior_set))
    ph_VIC_minus_Kent_and_EU = list(set(ph_VIC_prior_set) - set(ph_Kent_prior_set) & set(ph_EU_prior_set))
    h_sapiens = HSapiensManager(r"C:\studies\code\NetProp\data\H_sapiens.net")
    human_ppi = h_sapiens.get_data()
    problems = []
    for gene_id in ph_VIC_minus_Kent:
        try:
            human_ppi.nodes[gene_id][PropagationGraph.CONTAINER_PROPERTIES_KEY].source_of = {"liquid"}
        except KeyError:
            print(f"gene {gene_id} for prior set does not, in fact, appear to be a human gene")
            problems.append(gene_id)
    print(f"ignored the following genes: {problems}")
    ph_VIC_minus_Kent = [gene_id for gene_id in ph_VIC_minus_Kent if gene_id not in problems]
    problems = []
    for gene_id in ph_VIC_minus_Kent_and_EU:
        try:
            human_ppi.nodes[gene_id][PropagationGraph.CONTAINER_PROPERTIES_KEY].source_of = {"liquid"}
        except KeyError:
            print(f"gene {gene_id} for prior set does not, in fact, appear to be a human gene")
            problems.append(gene_id)
    print(f"ignored the following genes: {problems}")
    ph_VIC_minus_Kent_and_EU = [gene_id for gene_id in ph_VIC_minus_Kent_and_EU if gene_id not in problems]

    p1 = Process(target=generic_propagate_on_random_networks, args=(VIC_minus_Kent, r"C:\studies\code\NetProp\results\variants_blitz\super_randomized_networks", r"C:\studies\code\NetProp\results\variants_blitz\randomized_networks_propagations\VIC_minus_Kent", "VIC_minus_Kent"))
    p2 = Process(target=generic_propagate_on_random_networks, args=(VIC_minus_Kent_and_EU, r"C:\studies\code\NetProp\results\variants_blitz\super_randomized_networks", r"C:\studies\code\NetProp\results\variants_blitz\randomized_networks_propagations\VIC_minus_Kent_and_EU", "VIC_minus_Kent_and_EU"))
    p3 = Process(target=generic_propagate_on_random_networks, args=(ph_VIC_minus_Kent, r"C:\studies\code\NetProp\results\variants_blitz\super_randomized_networks", r"C:\studies\code\NetProp\results\variants_blitz\randomized_networks_propagations\ph_VIC_minus_Kent", "ph_VIC_minus_Kent"))
    p4 = Process(target=generic_propagate_on_random_networks, args=(ph_VIC_minus_Kent_and_EU, r"C:\studies\code\NetProp\results\variants_blitz\super_randomized_networks", r"C:\studies\code\NetProp\results\variants_blitz\randomized_networks_propagations\ph_VIC_minus_Kent_and_EU", "ph_VIC_minus_Kent_and_EU"))
    p1.start()
    p2.start()
    p3.start()
    p4.start()
    p1.join()
    p2.join()
    p3.join()
    p4.join()
