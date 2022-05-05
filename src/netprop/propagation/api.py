import json
from typing import Union
from multiprocessing import Process, Queue, cpu_count

from pydantic import BaseModel, Field

from netprop.networks.loaders.base import BaseNetworkLoader
from .api_helpers import _prop_from_conf, _prop_from_args_worker, _init_normalizing_propagaters
from netprop.models.config_models import ConfigModel, SuppressedSetModel
from netprop.models import ArgsModel


# This wrapper function is here to make it easy to integrate configuration duplicates in the future.
def propagate_from_config(config_path, ordering=None, max_processes=cpu_count() - 2):
    #TODO appears to be an internal bug in pydantic that doesn't allow it to read suppressed nodes - restore to commented line when it is resolved
    # config = ConfigModel.parse_file(config_path)
    with open(config_path, 'r') as handler:
        x = json.load(handler)
    config = ConfigModel.parse_obj(x)
    config.suppressed_set = [SuppressedSetModel.parse_obj(y) for y in x.get("suppressed_set", [])]
    if config.suppressed_set is None or len(config.suppressed_set) == 0:
        config.suppressed_set = SuppressedSetModel(id="no_knockouts")

    if not config.output_dir_path:
        config.output_dir_path = ""

    _prop_from_conf(config, ordering=ordering, max_processes=max_processes)


def prop_from_args(inputs: list[ArgsModel], loader_class: Union[str, BaseNetworkLoader], output_dir, max_processes=max(cpu_count() - 2, 1), **kwargs):
    queue = Queue()
    inputs = sorted(inputs, key=lambda inp: inp.network_init_args[0]) # Sort by network loader class so that workers reload networks infrequently
    for inp in inputs:
        queue.put((inp.prop_id, inp.network_init_args, inp.prior_set, inp.suppressed_set, inp.norm, inp.norm_kwargs))
    for i in range(max_processes):
        queue.put("HALT")
    processes = [Process(target=_prop_from_args_worker, args=(loader_class, output_dir, queue), kwargs=kwargs)
                         for _ in range(max_processes)]
    for p in processes:
        p.start()
    for p in processes:
        p.join()
