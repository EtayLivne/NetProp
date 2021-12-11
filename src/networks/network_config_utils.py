import re
import json
from typing import List
from pathlib import Path

from generic_utils.utils import listify
from networks.loaders.base import MultiNetworkLoader
from models.config_models import NetworksParametersModel, NetworkInclusionParametersModel


global_network_default_id_counter = 0


# TODO: better mechanism than the multi parameter, which is cumbersone and invites configuration errors
def _load(config: NetworksParametersModel, loader_classes):
    loader_class = loader_classes.get(config.loader_class, None)
    if not loader_class:
        raise ValueError(f"No network loader named {config.loader_class} detected")
    if config.multi != issubclass(loader_class, MultiNetworkLoader):
        raise ValueError(f"network {config.id}  multi value {config.multi} doesn't match network loader {config.loader_class}")
    loader = loader_class(*config.init_args, **config.init_kwargs)
    return loader.load()


def _set_metadata(network, config: NetworksParametersModel):
    if config.id:
        network.graph[config.id_key] = config.id
    else:
        network.graph[config.id_key] = global_network_default_id_counter
        global_network_default_id_counter += 1

    for k, v in config.metadata.items():
        network.graph[k] = v


def network_from_config(network_config, loader_classes):
    network = _load(network_config, loader_classes)
    _set_metadata(network, network_config)
    return network


def _unpack_multinetwork_config(network_config, loader_classes):
    return _load(network_config, loader_classes)


def _get_all_config_files(sources: List[str]):
    network_config_files = []
    for source in sources:
        path = Path(source)
        if path.is_dir():
            network_config_files.extend([f for f in path.glob("*.json")])
        elif path.is_file():
            network_config_files.append(path)
        else:
            raise ValueError(f"path {path} was given as a network config path, but does not exist")
    return network_config_files


def _file_passes_filter(config: NetworkInclusionParametersModel, file_path: str):
    if config.exclude_files:
        if any([re.search(re.escape(exclusion), str(file_path)) or exclusion == str(file_path) for
                exclusion in config.exclude_files]):
            return False
    elif config.include_files:
        if not any([re.search(re.escape(inclusion), str(file_path)) or inclusion == str(file_path) for
                    inclusion in config.include_files]):
            return False

    return True


def _id_passes_filter(config: NetworkInclusionParametersModel, network_id: str):
    if network_id:
        if config.exclude_ids and any([exclusion == network_id for exclusion in config.exclude_ids]):
            return False
        elif config.include_ids and not any([inclusion == network_id for inclusion in config.exclude_ids]):
            return False

    return True


def _add_to_unqiue_id_list(conf, unique_ids: dict, lst: list):
    if conf.id is not None and conf.id in unique_ids:
        raise ValueError(f"network cons with duplicate id {conf.id}")
    unique_ids[conf.id] = True
    lst.append(conf)


def single_network_config_iter(config: NetworkInclusionParametersModel, loader_classes):
    sources = config.path if isinstance(config.path, list) else [config.path]
    network_config_files = _get_all_config_files(sources)

    unpacked_configs = []
    unique_ids = {}
    for config_file in network_config_files:
        if not _file_passes_filter(config, config_file):
            continue

        with open(config_file, 'r') as handler:
            configs = listify(json.load(handler))

        # Handle multinetwork configs and create a flat list of single network configs, where each ID is unique
        for network_conf in configs:
            c = NetworksParametersModel.parse_obj(network_conf)
            if not _id_passes_filter(config, c.id):
                continue
            new_configs = _unpack_multinetwork_config(c, loader_classes) if c.multi else [c]

            for new_conf in new_configs:
                _add_to_unqiue_id_list(new_conf, unique_ids, unpacked_configs)
                if not [x for x in unpacked_configs if x.id == "combined_human_merged_covid"]:
                    x = 7

    for network_conf in unpacked_configs:
        yield network_conf
