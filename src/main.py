import json
import inspect
import sys
from os.path import join as os_path_join
from random import choice
from datetime import datetime
from importlib import import_module
from data_extractors import *
import common.network_loaders as network_loaders
from network_loader import BaseNetworkLoader
from propagater import Propagater, PropagationNetwork
from common.data_handlers.managers import HSapiensManager
from common.data_handlers.extractors import GeneInfoExtractor
from config_models import ConfigModel, PPIModel, PropagationSourceModel, PropagationTargetModel
from results_models import PropagationNetworkModel, PropagationResultModel, HaltConditionOptionModel,\
                           PropagationProteinModel
from constants import NodeAttrs


#TODO: find reasonable way to enforce set behavior on pyantics models. Currently cannot support sets of models so everything is a list
# reason for this is that .dict() will convert any member of the set to a dict, but a dict isn't hashable, which then crushes the program



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
        prop.deprecated_propagate(prior_set)
        node_tuples = sorted([(node, data[prop._LIQUIDS]["liquid"]) for node, data in prop._network.nodes(data=True)],
                             key=lambda x: x[1], reverse=True)
        output_path = output_dir + f"{output_file_name_iterator}"
        with open(output_path, 'w') as handler:
            json.dump([{"gene_id": node_tuple[0], "liquid": node_tuple[1], "in_prior_set": node_tuple[0] in prior_set} for node_tuple in node_tuples], handler, indent=4)
        print("yo")


def network_from_config(network_config, loader_classes):
    loader_class = loader_classes.get(network_config.loader_class, None)
    if not loader_class:
        raise ValueError(f"No network loader named {network_config.loader_class} detected")
    loader = loader_class(network_config.source_file)
    inputs = network_config.loader_input_dict or dict()

    propagation_network = loader.load_network(**inputs)
    propagation_network.graph[PPIModel.source_file_key] = network_config.source_file
    propagation_network.graph[PPIModel.id_key] = network_config.id
    propagation_network.graph[PPIModel.protein_id_class_key] = network_config.protein_id_class
    propagation_network.graph[PPIModel.directed_key] = network_config.directed
    return propagation_network


#TODO add support for suppressed nodes and variable halt conditions
def record_propagation_result(propagater, suppressed_nodes, file_path, propagation_id):
    network = PropagationNetworkModel(
        source_file=propagater.network.graph[PPIModel.source_file_key],
        network_id=propagater.network.graph[PPIModel.id_key],
        directed=propagater.network.graph[PPIModel.directed_key],
        suppressed_nodes=suppressed_nodes
    )

    max_steps = propagater.max_steps if propagater.max_steps != propagater.NO_MAX_STEPS else None
    min_gap = propagater.min_gap if propagater.min_gap != propagater.NO_MIN_GAP else None
    halt_conditions = HaltConditionOptionModel(max_steps=max_steps, min_gap=min_gap)

    nodes = [PropagationProteinModel(
        id=node,
        source_of=data[PropagationNetwork.CONTAINER_KEY].source_of,
        target_of=data[PropagationNetwork.CONTAINER_KEY].target_of,
        liquids=propagater.node_liquids(node),
        id_type=propagater.network.graph[PPIModel.protein_id_class_key],
        species=data[NodeAttrs.SPECIES_ID.value]
    ) for node, data in propagater.network.nodes(data=True) if node not in suppressed_nodes]

    propagation_result = PropagationResultModel(
        id=propagation_id,
        network=network,
        prior_set_confidence=propagater.source_confidence,
        halt_conditions=halt_conditions,
        nodes=nodes
    )

    with open(file_path, 'w') as json_handler:
        json_handler.write(propagation_result.json(indent=4, exclude_none=True))

# TODO add support for configurable halt condition
def propagate(network, propagation_config, output_dir):
    max_steps = propagation_config.halt_conditions.max_steps or Propagater.NO_MAX_STEPS
    min_gap = propagation_config.halt_conditions.min_gap or Propagater.NO_MIN_GAP

    prior_set_confidence = propagation_config.prior_set_confidence
    prior_set, target_set, suppressed_set\
        = propagation_config.prior_set, propagation_config.target_set, propagation_config.suppressed_set

    for prior, target in zip(prior_set, propagation_config.target_set):
        network.nodes[prior.id][PropagationNetwork.CONTAINER_KEY].source_of = prior.source_of
        network.nodes[target.id][PropagationNetwork.CONTAINER_KEY].target_of = target.target_of

    propagater = Propagater(network, prior_set_confidence, max_steps=max_steps, min_gap=min_gap)
    propagater.propagate([p.id for p in prior_set], suppressed_set=suppressed_set)

    propagation_id = propagation_config.id or network[PPIModel.network_id_key] + "_" + str(datetime.now())
    record_propagation_result(propagater,
                              suppressed_set,
                              os_path_join(output_dir, f"{propagation_id}.json"),
                              propagation_id)

    for prior, target in zip(prior_set, propagation_config.target_set):
        network.nodes[prior.id][PropagationNetwork.CONTAINER_KEY].source_of = set()
        network.nodes[target.id][PropagationNetwork.CONTAINER_KEY].target_of = set()


# def analytic_propagation_result()


def main():
    config_path = "example_config.json"
    config = ConfigModel.parse_file(config_path)
    _name, _cls = 0, 1
    loader_classes = {item[_name]: item[_cls] for item in inspect.getmembers(network_loaders, inspect.isclass)
                      if issubclass(item[_cls], BaseNetworkLoader)}
    network = network_from_config(config.ppi_config, loader_classes)
    output_dir = config.output_dir_path
    for propagation_config in config.propagations:
        propagate(network, propagation_config, output_dir)



if __name__ == "__main__":
    main()