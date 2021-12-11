import json
from functools import partial
from typing import Callable

from models.config_models import *


def save_as_json(model: BaseModel, path: str):
    with open(path, 'w') as handler:
        json.dump(model.dict(exclude_unset=True, exclude_defaults=True), handler, indent=4)


def copy_and_modify_config(original_config_path: str, mod_method: Callable, modified_conf_path: str):
    conf = ConfigModel.parse_file(original_config_path)
    mod_method(conf)
    save_as_json(conf, modified_conf_path)


def generic_modify(original_config: ConfigModel,
                   networks: NetworkInclusionParametersModel = None,
                   prior_set: list[PriorSetModel] = None,
                   suppressed_set: list[SuppressedSetModel] = None,
                   halt_conditions: HaltConditionOptionModel = None,
                   method: str = None,
                   prop_id: str = None,
                   output_dir_path: str = None):
    if networks is not None:
        original_config.networks = networks
    if prior_set is not None:
        original_config.prior_set = prior_set
    if suppressed_set is not None:
        original_config.suppressed_set = suppressed_set
    if halt_conditions is not None:
        original_config.halt_conditions = halt_conditions
    if method is not None:
        original_config.method = method
    if prop_id is not None:
        original_config.id = prop_id
    if output_dir_path is not None:
        original_config.output_dir_path = output_dir_path




# import networkx as nx
# from networks.loaders.single import CombinedHumanCovidNetworkLoader
# network = CombinedHumanCovidNetworkLoader(r"C:\studies\thesis\code\NetProp\data\H_sapiens\H_spaiens_nov_2021.net",
#                                           r"C:\studies\thesis\code\NetProp\data\ndex_covid_human_ppi.cx",
#                                           r"C:\studies\thesis\code\NetProp\data\symbol_to_entrezgene_2021.json").load()
#
# knockout_candidates = []
# for n in nx.ego_graph(network, "covid", center=False, radius=1).nodes():
#     knockout_candidates.append(SuppressedSetModel(id=f"{n}_knockout", nodes=[n]))
#
#
# original_conf_path = r"C:\studies\thesis\code\temp\configurations\new\example_config.json"
# output_conf_path = r"C:\studies\thesis\code\temp\configurations\new\covid_neighbors_knockout_conf.json"
#
# mod = partial(generic_modify, suppressed_set=knockout_candidates)
#
# copy_and_modify_config(original_conf_path, mod, output_conf_path)