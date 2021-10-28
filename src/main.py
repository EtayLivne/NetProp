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
from config_models import ConfigModel, PPIModel, PropagationParametersModel
from results_models import PropagationNetworkModel, PropagationResultModel, HaltConditionOptionModel,\
                           PropagationProteinModel
from constants import NodeAttrs
import sys


#TODO: find reasonable way to enforce set behavior on pyantics models. Currently cannot support sets of models so everything is a list
# reason for this is that .dict() will convert any member of the set to a dict, but a dict isn't hashable, which then crushes the program


def generate_rank_equivalent_random_network(network):
    switches = 100 * len(list(network.edges))
    switch_types = ["place", "target", "source"]
    edge_source_sample = []
    for node in network.nodes:
        edge_source_sample.extend([node] * network.degree(node))

    while switches:
        switch_type = choice(switch_types)
        source_1 = choice(edge_source_sample)
        target_1 = choice(list(network.neighbors(source_1)))
        edge_1 = (source_1, target_1)
        edge_1_attrs = dict(network.edges[edge_1])

        # given the first random edge, generate a second one where the change would be an actual change
        while True:
            source_2 = choice(edge_source_sample)
            target_2 = choice(list(network.neighbors(source_2)))
            if (switch_type == "place" and (source_2 not in edge_1)) or\
               (switch_type == "source" and (source_2 not in edge_1 and source_1 != target_2)) or\
               (switch_type == "target" and (target_2 not in edge_1 and source_2 != target_1)):
                break
        edge_2 = (source_2, target_2)
        edge_2_attrs = dict(network.edges[edge_2])

        network.remove_edge(*edge_1)
        network.remove_edge(*edge_2)
        if switch_type == "place":
            network.add_edge(*edge_1, **edge_2_attrs)
            network.add_edge(*edge_2, **edge_1_attrs)
        if switch_type == "target":
            network.add_edge(edge_1[0], edge_2[1], **edge_1_attrs)
            network.add_edge(edge_2[0], edge_1[1], **edge_2_attrs)
        if switch_type == "source":
            network.add_edge(edge_2[0], edge_1[1], **edge_1_attrs)
            network.add_edge(edge_1[0], edge_2[1], **edge_2_attrs)

        switches -= 1


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

    prior_set_confidence, method = propagation_config.prior_set_confidence, propagation_config.method
    prior_set, target_set, suppressed_set\
        = propagation_config.prior_set, propagation_config.target_set, propagation_config.suppressed_set

    for prior in prior_set:
        network.nodes[prior.id][PropagationNetwork.CONTAINER_KEY].source_of = prior.source_of

    for target in propagation_config.target_set:
        network.nodes[target.id][PropagationNetwork.CONTAINER_KEY].target_of = target.target_of

    propagater = Propagater(network, prior_set_confidence, max_steps=max_steps, min_gap=min_gap)
    propagater.propagate([p.id for p in prior_set], suppressed_set=suppressed_set, propagation_method=method)

    propagation_id = propagation_config.id or network[PPIModel.network_id_key] + "_" + str(datetime.now())
    record_propagation_result(propagater,
                              suppressed_set,
                              os_path_join(output_dir, f"{propagation_id}.json"),
                              propagation_id)

    for prior, target in zip(prior_set, propagation_config.target_set):
        network.nodes[prior.id][PropagationNetwork.CONTAINER_KEY].source_of = set()
        network.nodes[target.id][PropagationNetwork.CONTAINER_KEY].target_of = set()


# def analytic_propagation_result()

def load_config(config_path):
    config = ConfigModel.parse_file(config_path)
    if config.global_propagation_params:
        global_params = config.global_propagation_params.dict(exclude_unset=True)
        propagation_dicts = [merge_dicts(global_params, p.dict(exclude_unset=True)) for p in config.propagations]
        config.propagations = [PropagationParametersModel(**d) for d in propagation_dicts]

    if val_error := next(iter([p for p in config.propagations if not p.validate_completeness()]), None):
        error_propagation_id = val_error.id or str(val_error)
        raise ValueError(f"propagation identified as {error_propagation_id} doesnt have every required field")

    return config

def merge_dicts(*dicts):
    new_dict = dict()
    for d in dicts:
        new_dict.update(d)
    return new_dict

def main():
    config_path = sys.argv[1]
    config = load_config(config_path)

    _name, _cls = 0, 1
    loader_classes = {item[_name]: item[_cls] for item in inspect.getmembers(network_loaders, inspect.isclass)
                      if issubclass(item[_cls], BaseNetworkLoader)}
    network = network_from_config(config.ppi_config, loader_classes)
    output_dir = config.output_dir_path
    for propagation_config in config.propagations:
        propagate(network, propagation_config, output_dir)


if __name__ == "__main__":
    main()

    # network = network_loaders.HumanCovidHybridNetworkLoader("C:\\studies\\code\\NetProp\\data\\H_sapiens.net").load_network()
    # prior_set = [
    #     {
    #         "id": node,
    #         "source_of": [node]
    #     } for node, data in network.nodes(data=True) if data["species_id"] == "sars-cov-2"
    # ]
    #
    # with open("temp\\conf.json", 'w') as handler:
    #     json.dump({
    #         "ppi_config": {
    #             "id": "mr_dummy",
    #             "source_file": "C:\\studies\\code\\NetProp\\data\\H_sapiens.net",
    #             "loader_class": "HumanCovidHybridNetworkLoader",
    #             "protein_id_class": "entrezgene"
    #         },
    #         "global_propagation_params": {
    #             "prior_set_confidence": 0.8,
    #             "halt_conditions": {
    #                 "min_gap": 1e-4
    #             }
    #         },
    #         "propagations": [
    #             {
    #                 "id": "all_covid_sources",
    #                 "prior_set": prior_set
    #             }
    #         ],
    #         "output_dir": "temp\\propagation_results"
    #     }, handler, indent=4)