import inspect
import sys
from os.path import join as os_path_join
from random import choice
from multiprocessing import Queue, Pool
from datetime import datetime
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
from os import path as os_path

#TODO: find reasonable way to enforce set behavior on pyantics models. Currently cannot support sets of models so everything is a list
# reason for this is that .dict() will convert any member of the set to a dict, but a dict isn't hashable, which then crushes the program


def generate_rank_equivalent_random_network(network, edge_switch_factor):
    switches = edge_switch_factor * len(list(network.edges))
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
            edge_2 = (source_2, target_2)
            if (switch_type == "place" and (len(set(edge_1 + edge_2)) > 2)) or\
               (switch_type != "place" and (len(set(edge_1 + edge_2)) == 4)):
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


def worker_main(ppi_config, output_dir, propagation_queue):
    _name, _cls = 0, 1
    loader_classes = {item[_name]: item[_cls] for item in inspect.getmembers(network_loaders, inspect.isclass)
                      if issubclass(item[_cls], BaseNetworkLoader)}
    network = network_from_config(ppi_config, loader_classes)

    while True:
        propagation_config = propagation_queue.get(block=True)
        if propagation_config == "HALT":
            break
        print(f"Now propagating {propagation_config.id}")
        propagate(network, propagation_config, output_dir)


def randomizer_main(original_network_file, network_loader_class, queue):
    network_loader = network_loader_class(original_network_file)


    while True:
        randomization_request = queue.get(block=True)

        if randomization_request == "HALT":
            break
        output_path = randomization_request["output_path"]
        switch_factor = randomization_request["edge_switch_factor"]
        print(f"now handling request for {output_path}")
        network = network_loader.load_network()
        generate_rank_equivalent_random_network(network, switch_factor)
        network_loader_class.record_network(network, output_path)


def main():
    config_path = sys.argv[1]
    config = load_config(config_path)
    propagations_queue = Queue()
    worker_pool = Pool(10, worker_main, (config.ppi_config, config.output_dir_path, propagations_queue))
    for propagations_config in config.propagations:
        propagations_queue.put(propagations_config)
    for worker in range(10):
        propagations_queue.put("HALT")

    propagations_queue.close()
    propagations_queue.join_thread()
    worker_pool.close()
    worker_pool.join()


def randomize_network(number_of_randomizations, edge_switch_factor, original_network_file,
                      network_loader_class, output_dir):
    queue = Queue()
    randomizers = Pool(10, randomizer_main, (original_network_file, network_loader_class, queue))
    original_network_name = os_path.basename(original_network_file).split(".")[0]
    for i in range(number_of_randomizations):
        queue.put({"edge_switch_factor": edge_switch_factor,
                   "output_path": os_path.join(output_dir, original_network_name + f"_{i}")})

    for i in range(10):
        queue.put("HALT")

    queue.close()
    queue.join_thread()
    randomizers.close()
    randomizers.join()

if __name__ == "__main__":
    main()
    # randomize_network(100, 50, r"C:\studies\code\NetProp\data\H_sapiens.net", network_loaders.HSapeinsNetworkLoader,
    #                   r"C:\studies\code\NetProp\data\randomized_h_sapiens_with_covid")



    # with open(r"C:\studies\code\NetProp\src\temp\configurations\original_h_sapiens_covid_conf.json", "r") as handler:
    #     conf = json.load(handler)
    # conf["output_dir_path"] = r"C:\studies\code\NetProp\src\temp\propagations\gene_knockouts"
    # propagations = []
    # full_network_propagation = conf["propagations"][0]
    # network = network_loaders.HumanCovidHybridNetworkLoader(conf["ppi_config"]["source_file"]).load_network()
    # for n in [n for n in network.nodes if network.nodes[n]["species_id"] == "human"][1000:4000]:
    #     knockout_propagation = full_network_propagation.copy()
    #     knockout_propagation.update({"suppressed_set": [n], "id": f"{n}_knockout"})
    #     propagations.append(knockout_propagation)
    #
    # conf["propagations"] = propagations
    # print(f"{len(propagations)} single knockouts")
    # with open("temp/configurations/single_knockouts_conf.json", "w") as handler:
    #     json.dump(conf, handler, indent=4)


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
