
from random import choice
from importlib import import_module
from data_extractors import *
from propagater import Propagater
from common.data_handlers.managers import HSapiensManager
from common.data_handlers.extractors import GeneInfoExtractor
from config_models import ConfigModel

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


def network_form_config(network_config):
    loader_class = getattr(import_module(network_config.network_loader), "NetworkLoader")
    inputs = network_config.network_loader_input_dict or dict()
    propagation_network = loader_class.load_network(network_config.source_file, **inputs)




def main():
    config_path = "config.json"
    config = ConfigModel.parse_file(config_path)
    network = network_from_config(config.ppi_config)
    output_dir = config.output_dir_path
    for propagation in config.propagations:
        propagate(output_dir, propagation)


if __name__ == "__main__":
    main()