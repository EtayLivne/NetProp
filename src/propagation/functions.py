import inspect
from pathlib import Path
from datetime import datetime
from multiprocessing import Pool, Queue, cpu_count

from utils.utils import listify
from networks.network_loader import BaseNetworkLoader
import networks.single_network_loaders as network_loaders
import networks.multi_network_loaders as multi_network_loaders
from propagation.classes import Propagater, PropagationNetwork
from models.config_models import ConfigModel, NetworksParametersModel, SuppressedSetModel
from networks.network_config_utils import single_network_config_iter, network_from_config
from models.results_models import PropagationResultModel, HaltConditionOptionModel, PropagationNodeModel


def record_propagation_result(propagater, suppressed_nodes, file_path, propagation_id):
    network = propagater.network.graph
    max_steps = propagater.max_steps if propagater.max_steps != propagater.NO_MAX_STEPS else None
    min_gap = propagater.min_gap if propagater.min_gap != propagater.NO_MIN_GAP else None
    halt_conditions = HaltConditionOptionModel(max_steps=max_steps, min_gap=min_gap)
    non_metadata_fields = Propagater.propagation_related_node_metadata()
    nodes = {node: PropagationNodeModel(id=node,
                                        source_of=data[PropagationNetwork.CONTAINER_KEY].source_of,
                                        liquids=propagater.node_liquids(node),
                                        metadata={k: v for k, v in data.items() if k not in non_metadata_fields})
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


def propagate(network, config: ConfigModel, output_dir):
    max_steps = config.halt_conditions.max_steps or Propagater.NO_MAX_STEPS
    min_gap = config.halt_conditions.min_gap or Propagater.NO_MIN_GAP

    method = config.method
    suppressed_set = config.suppressed_set.nodes
    prior_set, prior_set_confidence = config.prior_set.nodes, config.prior_set.confidence
    for prior in prior_set:
        network.nodes[prior.id][PropagationNetwork.CONTAINER_KEY].source_of = prior.source_of

    network.graph["prior_set"] = prior_set
    network.graph["suppressed_nodes"] = suppressed_set

    propagater = Propagater(network, prior_set_confidence, max_steps=max_steps, min_gap=min_gap)
    propagater.propagate([p.id for p in prior_set], suppressed_set=suppressed_set, propagation_method=method)

    propagation_id = config.id or network.graph[NetworksParametersModel.id_key] + "_" + str(datetime.now())
    record_propagation_result(propagater,
                              suppressed_set,
                              str(Path(output_dir) / f"{propagation_id}.json"),
                              propagation_id)

    for prior in prior_set:
        network.nodes[prior.id][PropagationNetwork.CONTAINER_KEY].source_of = set()


def propagation_worker(network_conf, loader_classes, queue):
    network = network_from_config(network_conf, loader_classes)

    while True:
        propagation_config = queue.get(block=True)
        if propagation_config == "HALT":
            break
        print(f"Now propagating {propagation_config.id}")
        propagate(network, propagation_config, propagation_config.output_dir_path)


def launch_multiprocess_propagation(config: ConfigModel, max_processes=cpu_count()):
    _name, _cls = 0, 1
    loader_classes = {item[_name]: item[_cls] for item in
                      inspect.getmembers(network_loaders, inspect.isclass) + inspect.getmembers(multi_network_loaders, inspect.isclass)
                      if issubclass(item[_cls], BaseNetworkLoader)}

    suppressed_nodes_sets = listify(config.suppressed_set)
    prior_sets = listify(config.prior_set)
    network_replicates = len(suppressed_nodes_sets) * len(prior_sets)
    pool_size = min(network_replicates, max_processes)
    network_counter = 0
    for network_conf in single_network_config_iter(config.networks, loader_classes):
        if not network_conf.id:
            network_conf.id = f"network_{network_counter}"
            network_counter += 1
        network_dir = Path(config.output_dir_path) / str(network_conf.id)
        network_dir.mkdir(exist_ok=True)
        queue = Queue()
        specific_config = config.copy()
        prior_set_counter = 0
        suppressed_nodes_counter = 0
        for prior_set in prior_sets:
            if not prior_set.id:
                prior_set.id = f"prior_set_{prior_set_counter}"
                prior_set_counter += 1
            prior_set_dir = network_dir / str(prior_set.id)
            prior_set_dir.mkdir(exist_ok=True)

            for suppressed_nodes in suppressed_nodes_sets:
                if not suppressed_nodes.id:
                    suppressed_nodes.id = f"knockout_set_{suppressed_nodes_counter}"
                    suppressed_nodes_counter += 1

                specific_config.prior_set = prior_set
                specific_config.output_dir_path = prior_set_dir
                specific_config.suppressed_set = suppressed_nodes
                specific_config.id = suppressed_nodes.id

                queue.put(specific_config)

        for i in range(pool_size):
            queue.put("HALT")

        worker_pool = Pool(pool_size, propagation_worker, (network_conf, loader_classes, queue))
        queue.close()
        queue.join_thread()
        worker_pool.close()
        worker_pool.join()


# This wrapper function is here to make it easy to integrate configuration duplicates in the future.
def propagate_from_config(config_path):
    config = ConfigModel.parse_file(config_path)
    if config.suppressed_set is None:
        config.suppressed_set = SuppressedSetModel(id="no_knockouts")
    if not config.output_dir_path:
        config.output_dir_path = Path.cwd()

    launch_multiprocess_propagation(config)


if __name__ == "__main__":
    propagate_from_config(r"/temp/configurations/new/dummy.json")