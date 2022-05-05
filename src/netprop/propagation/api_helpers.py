from typing import Union, Callable
from multiprocessing import Queue, Pool, cpu_count, Process
from pathlib import Path
import networkx as nx
from datetime import datetime

from netprop.models import ConfigModel, HaltConditionOptionModel, PropagationResultModel, PropagationNodeModel, NetworksParametersModel
from netprop.networks.loaders.base import BaseNetworkLoader, loader_registry
from netprop.networks.network_config_utils import single_network_config_iter, network_from_config, listify


from .classes.propagater import Propagater
from .classes.normalizing_propagater import norm_prop_registry

NETWORK_ORDERING_KEYWORD = "network"
PRIOR_SET_ORDERING_KEYWORD = "prior_set"
SUPPRESSED_NODES_ORDERING_KEYWORD = "suppressed_nodes"

def _init_normalizing_propagaters():
    from netprop.propagation.classes import rdpn_normalizing_propagater


def _prop_from_conf(config: ConfigModel, max_processes=cpu_count() - 2, ordering=None):
    if ordering is None:
        ordering = {NETWORK_ORDERING_KEYWORD: 0, PRIOR_SET_ORDERING_KEYWORD: 1, SUPPRESSED_NODES_ORDERING_KEYWORD: 2}

    _name, _cls = 0, 1

    suppressed_nodes_sets = listify(config.suppressed_set)
    prior_sets = listify(config.prior_set)
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

                output_path = _prop_from_conf_determine_output_path(ordering,
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
            for i in range(max_processes):
                queue.put("HALT")

            _consume_task_queue(_prop_from_conf_worker, queue, max_processes)
            queue_size = 0
            queue = Queue()

    if queue_size:
        for i in range(max_processes):
            queue.put("HALT")

        _consume_task_queue(_prop_from_conf_worker, queue, max_processes)


def _consume_task_queue(task_consuming_func: Callable, task_queue: Queue, max_processes: int) -> None:
    if max_processes > 1:
        worker_pool = [Process(target=task_consuming_func, args=(task_queue,)) for i in range(max_processes)]
        for p in worker_pool:
            p.start()
        for p in worker_pool:
            p.join()
    else:
        task_consuming_func(task_queue)
def _prop_from_conf_determine_output_path(ordering: dict,
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


def _prop_from_conf_worker(queue):
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
        _start_prop_from_conf(network, propagation_config, propagation_config.output_dir_path)


def _start_prop_from_conf(network, config: ConfigModel, output_dir):
    if config.suppressed_set:
        sup_set = config.suppressed_set.nodes
    else:
        sup_set = None
    _start_prop_from_args(network, config.method, sup_set, {p.id: p.source_of for p in config.prior_set.nodes},
                          config.prior_set.confidence, config.id, config.norm, config.norm_kwargs, output_dir,
                          max_steps=config.halt_conditions.max_steps, min_gap=config.halt_conditions.min_gap)


def _prop_from_args_worker(nw_load_class: Union[str, BaseNetworkLoader], output_dir: str, queue: Queue, **kwargs):
    method = kwargs.get("method", "iterative")
    prior_set_confidence = kwargs.get("prior_set_confidence", 0.7)
    min_gap = kwargs.get("min_gap", None)
    max_steps = kwargs.get("max_steps", 20)

    curr_nw_args = None
    curr_nw = None
    while True:
        task = queue.get()
        if task == "HALT":
            break
        else:
            prop_id, nw_args, prior_set, suppressed_set, norm, norm_kwargs = task
            if nw_args != curr_nw_args:
                loader = nw_load_class if issubclass(nw_load_class, BaseNetworkLoader) else loader_registry[nw_load_class]
                curr_nw = loader(*nw_args).load()
                curr_nw_args = nw_args

        _start_prop_from_args(curr_nw,
                              method, suppressed_set, prior_set, prior_set_confidence, prop_id, norm, norm_kwargs,
                              output_dir,
                              max_steps=max_steps, min_gap=min_gap)

    return


def _start_prop_from_args(network: nx.Graph,
                          method: str, suppressed_set: Union[list[str], set[str]],
                          prior_set: dict[str, str], prior_set_confidence,
                          prop_id: str, norm: str, norm_kwargs: dict, output_dir: str,
                          max_steps: int=None, min_gap: float=None) -> None:

    if type(prior_set) is dict:
        for prior, source_of in prior_set.items():
            network.nodes[prior]["source_of"] = source_of
        priors = list(prior_set.keys())
    else:
        priors = prior_set

    if norm is not None:
        from netprop.propagation.classes import rdpn_normalizing_propagater
        prop_class = norm_prop_registry[norm]
        prop = prop_class(network, prior_set_confidence, priors, suppressed_set,
                      min_gap=min_gap, max_steps=max_steps)
        prop.propagate(norm_kwargs=norm_kwargs, normalize=True)
    else:
        prop = Propagater(network, prior_set_confidence, priors, suppressed_set,
                          min_gap=min_gap, max_steps=max_steps)

        prop.propagate()

    propagation_id = prop_id or network.graph[NetworksParametersModel.id_key] + "_" + str(datetime.now())
    record_propagation_result(prop,
                              str(Path(output_dir) / f"{propagation_id}.json"),
                              propagation_id)


def record_propagation_result(propagater, file_path, propagation_id):
    network = propagater.network.graph
    halt_conditions = HaltConditionOptionModel(max_steps=propagater.max_steps, min_gap=propagater.min_gap)
    state = propagater.get_network_state()
    nodes = {node: PropagationNodeModel(id=node,
                                        source_of=data.get("source_of", set()),
                                        liquids=state.loc[node].to_dict(),
                                        metadata=data) # Also records the source of field again. Allowing now for speed
             for node, data in propagater.network.nodes(data=True)}

    propagation_result = PropagationResultModel(
        id=propagation_id,
        network=network,
        prior_set_confidence=propagater.source_confidence,
        halt_conditions=halt_conditions,
        nodes=nodes
    )

    with open(file_path, 'w') as json_handler:
        json_handler.write(propagation_result.json(indent=4, exclude_none=True))