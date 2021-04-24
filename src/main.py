from random import choice
from itertools import cycle
from multiprocessing import Process

from data_extractors import *
from common.data_handlers.managers import HSapiensManager
from common.data_handlers.extractors import GeneInfoExtractor
from propagater import Propagater, PropagationContainer, PropagationGraph


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
            if source_2 in random_edge_1 or target_2 in random_edge_1:
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


def generic_propagate_on_random_networks(prior_set, prior_set_confidence,
                                         randomized_networks, output_dir, output_file_name_iterator):

    for randomized_network in randomized_networks:
        prop = Propagater(randomized_network, prior_set_confidence)
        prop.propagate(prior_set)
        node_tuples = sorted([(node, data[prop._LIQUIDS]["liquid"]) for node, data in prop.graph.nodes(data=True)],
                             key=lambda x: x[1], reverse=True)
        output_path = output_dir + f"{output_file_name_iterator}"
        with open(output_path, 'w') as handler:
            json.dump([{"gene_id": node_tuple[0], "liquid": node_tuple[1], "in_prior_set": node_tuple[0] in prior_set} for node_tuple in node_tuples], handler, indent=4)
        print("yo")


def load_propagation_network(source_file_path):
    try:
        manager = HSapiensManager(source_file_path)
    except:
        raise (f"failed to read network from file {source_file_path}")

    return manager.get_data()

# TODO april make this the core of a generic way to run propagations!
def propagation_iter(propagaters_iter, prior_sets_iter, prior_set_confidences_iter=Propagater.DEFAULT_CONFIDENCE_COEF, knockout_sets_iter=None):
    stopped_iterators = {}
    knockout_sets_iter = knockout_sets_iter or cycle([None])
    while True:
        propagater = next(propagaters_iter, "stopped")
        knockout_set = next(knockout_sets_iter, "stopped")
        prior_set = next(prior_sets_iter, "stopped")
        prior_set_confidence = next(prior_set_confidences_iter, "stopped")
        stopped_iterators = {i for i in [knockout_set, prior_set, prior_set_confidence] if i == "stopped"}
        if stopped_iterators:
            break
        yield propagater.propagate(prior_set, suppressed_nodes=knockout_set, normalize_flow=False)

    print(f'stopped iteration because the following iterators ran out: {stopped_iterators}')


def propagate_on_network(network_file_path, prior_sets, prior_set_confidences, output_dir, knockout_sets=None):
    network = HSapiensManager(file_path=network_file_path).get_data()
    prop = Propagater(network, 0)
    for i in range(min(len(prior_sets)))
    results_manager = prop.propagate(prior_set, suppressed_nodes=knockout_set)
    results_manager.dump_to_file(output_file_path)





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
    #     output_file_path = output_dir + f"\\ph_EU120_{randomized_human_ppi_file}.json"
    #     with open(output_file_path, 'w') as handler:
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
    #     output_file_path = output_dir + f"\\ph_Kent_{randomized_human_ppi_file}.json"
    #     with open(output_file_path, 'w') as handler:
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
    #     output_file_path = output_dir + f"\\ph_VIC_{randomized_human_ppi_file}.json"
    #     with open(output_file_path, 'w') as handler:
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
