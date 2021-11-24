import inspect
from os.path import join as os_path_join
from datetime import datetime
from multiprocessing import Pool, Queue
from propagater import Propagater, PropagationNetwork

from config_models import ConfigModel, PPIModel, PropagationParametersModel
from results_models import PropagationNetworkModel, PropagationResultModel, HaltConditionOptionModel,\
                           PropagationProteinModel
from constants import NodeAttrs
from network_loader import BaseNetworkLoader
import common.network_loaders as network_loaders


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


# TODO add support for suppressed nodes and variable halt conditions
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

    nodes = {node: PropagationProteinModel(id=node,
                                           source_of=data[PropagationNetwork.CONTAINER_KEY].source_of,
                                           target_of=data[PropagationNetwork.CONTAINER_KEY].target_of,
                                           liquids=propagater.node_liquids(node),
                                           id_type=propagater.network.graph[PPIModel.protein_id_class_key],
                                           species=data[NodeAttrs.SPECIES_ID.value])
             for node, data in propagater.network.nodes(data=True) if node not in suppressed_nodes}

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


def propagater_process_main(ppi_config, propagation_queue):
    _name, _cls = 0, 1
    loader_classes = {item[_name]: item[_cls] for item in inspect.getmembers(network_loaders, inspect.isclass)
                      if issubclass(item[_cls], BaseNetworkLoader)}
    network = network_from_config(ppi_config, loader_classes)

    while True:
        propagation_config = propagation_queue.get(block=True)
        if propagation_config == "HALT":
            break
        print(f"Now propagating {propagation_config.id}")
        propagate(network, propagation_config, propagation_config.output_dir_path)


def launch_multiprocess_propgation(config: ConfigModel, num_processes=10):
    propagations_queue = Queue()
    worker_pool = Pool(num_processes, propagater_process_main, (config.ppi_config, propagations_queue))
    for propagations_config in config.propagations:
        propagations_queue.put(propagations_config)
    for worker in range(num_processes):
        propagations_queue.put("HALT")

    propagations_queue.close()
    propagations_queue.join_thread()
    worker_pool.close()
    worker_pool.join()


def multiprocessed_propagate_from_config(config_path):
    config = load_config(config_path)
    launch_multiprocess_propgation(config)