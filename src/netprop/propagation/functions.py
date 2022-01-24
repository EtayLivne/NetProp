import inspect
from pathlib import Path
from datetime import datetime
from multiprocessing import Pool, Queue, cpu_count

from netprop.generic_utils.utils import listify
from .classes import Propagater, PropagationNetwork
from netprop.networks.loaders.base import BaseNetworkLoader
import netprop.networks.loaders.single as network_loaders
import netprop.networks.loaders.multi as multi_network_loaders
from netprop.models.config_models import ConfigModel, NetworksParametersModel, SuppressedSetModel
from netprop.networks.network_config_utils import single_network_config_iter, network_from_config
from netprop.models.results_models import PropagationResultModel, HaltConditionOptionModel, PropagationNodeModel

NETWORK_ORDERING_KEYWORD = "network"
PRIOR_SET_ORDERING_KEYWORD = "prior_set"
SUPPRESSED_NODES_ORDERING_KEYWORD = "suppressed_nodes"


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
             for node, data in propagater.network.nodes(data=True)} # used to have a "if node not in suppressed_nodes" clause

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


def propagation_worker(queue):
    network = None
    network_id = None
    while True:
        task = queue.get(block=True)
        if task == "HALT":
            break

        propagation_config, network_conf = task
        if network_conf.id != network_id:
            network = network_from_config(network_conf)
            network_id = network_conf.id

        print(f"Now propagating {propagation_config.output_dir_path} - {propagation_config.id}")
        propagate(network, propagation_config, propagation_config.output_dir_path)


def determine_output_path(ordering: dict,
                          network_name: str, prior_set_name: str, suppressed_set_name: str,
                          root_output_path: str):

    ordered_names = sorted([(NETWORK_ORDERING_KEYWORD, network_name),
                            (PRIOR_SET_ORDERING_KEYWORD, prior_set_name),
                            (SUPPRESSED_NODES_ORDERING_KEYWORD, suppressed_set_name)],
                           key=lambda tup: ordering.get(tup[0], -1))
    filtered_ordered_names = [tup[1] for tup in ordered_names if ordering.get(tup[0], -1) > 0]
    output_dir = Path(root_output_path) / "/".join(filtered_ordered_names[:-1])
    Path.mkdir(output_dir, exist_ok=True, parents=True)
    return output_dir / filtered_ordered_names[-1]


def launch_multiprocess_propagation(config: ConfigModel, max_processes=cpu_count() - 2, ordering=None):
    if ordering is None:
        ordering = {NETWORK_ORDERING_KEYWORD: 0, PRIOR_SET_ORDERING_KEYWORD: 1, SUPPRESSED_NODES_ORDERING_KEYWORD: 2}

    _name, _cls = 0, 1

    suppressed_nodes_sets = listify(config.suppressed_set)
    prior_sets = listify(config.prior_set)
    pool_size = max_processes
    network_counter = 0
    queue_size = 0
    queue = Queue()
    for network_conf in single_network_config_iter(config.networks):
        if not network_conf.id:
            network_conf.id = f"network_{network_counter}"
            network_counter += 1

        prior_set_counter = 0
        suppressed_nodes_counter = 0
        for prior_set in prior_sets:
            if not prior_set.id:
                prior_set.id = f"prior_set_{prior_set_counter}"
                prior_set_counter += 1

            for suppressed_nodes in suppressed_nodes_sets:
                if suppressed_nodes and not suppressed_nodes.id:
                    suppressed_nodes.id = f"knockout_set_{suppressed_nodes_counter}"
                    suppressed_nodes_counter += 1

                output_path = determine_output_path(ordering,
                                                    network_conf.id, prior_set.id, suppressed_nodes.id,
                                                    config.output_dir_path)
                specific_config = config.copy()
                specific_config.prior_set = prior_set
                specific_config.suppressed_set = suppressed_nodes
                specific_config.output_dir_path = str(output_path.parent)
                specific_config.id = str(output_path.name)

                queue.put((specific_config, network_conf))
                queue_size += 1

        if queue_size > 1000:
            for i in range(pool_size):
                queue.put("HALT")

            worker_pool = Pool(pool_size, propagation_worker, (queue,))
            queue.close()
            queue.join_thread()
            worker_pool.close()
            worker_pool.join()
            queue_size = 0
            queue = Queue()

    if queue_size:
        for i in range(pool_size):
            queue.put("HALT")

        worker_pool = Pool(pool_size, propagation_worker, (queue,))
        queue.close()
        queue.join_thread()
        worker_pool.close()
        worker_pool.join()

# This wrapper function is here to make it easy to integrate configuration duplicates in the future.
def propagate_from_config(config_path, ordering=None):
    import json
    #TODO appears to be an internal bug in pydantic that doesn't allow it to read suppressed nodes - restore to commented line when it is resolved
    # config = ConfigModel.parse_file(config_path)
    with open(config_path, 'r') as handler:
        x = json.load(handler)
    config = ConfigModel.parse_obj(x)
    config.suppressed_set = [SuppressedSetModel.parse_obj(y) for y in x.get("suppressed_set", [])]
    if config.suppressed_set is None or len(config.suppressed_set) == 0:
        config.suppressed_set = SuppressedSetModel(id="no_knockouts")

    # attach path to volume root
    if not config.output_dir_path:
        config.output_dir_path = ""
    # else:
    #     config.output_dir_path = attach_to_root(config.output_dir_path)

    launch_multiprocess_propagation(config, ordering=ordering)