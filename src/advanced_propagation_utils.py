from pydantic import BaseModel
from typing import Optional, Dict, Set
from config_models import ConfigModel, PropagationParametersModel
from propagation_utils import launch_multiprocess_propgation, load_config
import os
from pathlib import Path

class MatrixConfigModel(BaseModel):
    singular_config_path: str
    randomized_networks_dir: Optional[str]
    suppressed_duplicates: Optional[Dict[str, Set[str]]]
    output_root_dir: str


def duplicate_configurations_with_suppressed_nodes(propagation_list, suppressed_sets: dict):

    for propagation in propagation_list:
        new_propagations = [propagation.copy(deep=True) for _ in suppressed_sets]
        suppressed_set_names = list(suppressed_sets.keys())
        for i in range(len(new_propagations)):
            new_propagation = new_propagations[i]
            suppressed_set_name = suppressed_set_names[i]
            new_propagation.suppressed_set.update(suppressed_sets[suppressed_set_name])
            new_propagation.id = new_propagation.id + "_" + suppressed_set_name
            # appends even if the new suppression didn't add any new suppressed nodes,
            # to create predictably structured output
    propagation_list.extend(new_propagations)

def propagate_a_lot(matrix_config_path: str):
    config = MatrixConfigModel.parse_file(matrix_config_path)
    singular_config = load_config(config.singular_config_path)
    suppressed_duplicates = config.suppressed_duplicates or dict()

    duplicate_configurations_with_suppressed_nodes(singular_config.propagations, suppressed_duplicates)
    for propagation in singular_config.propagations:
        propagation.output_dir_path = config.output_root_dir

    # propagate on original networks and suppressed networks
    launch_multiprocess_propgation(singular_config)

    if config.randomized_networks_dir:
        for propagation in singular_config.propagations:
            os.mkdir(os.path.join(config.output_root_dir, propagation.id + "_randomized"))

        current_id_to_original_id = {}

        for file in os.listdir(config.randomized_networks_dir):
            for i in range(len(singular_config.propagations)):
                propagation = singular_config.propagations[i]
                propagation.output_dir_path = os.path.join(config.output_root_dir,
                                                           propagation.id + "_randomized")
                current_id_to_original_id[i] = propagation.id
                propagation.id = file

            singular_config.ppi_config.source_file = os.path.join(config.randomized_networks_dir, file)
            singular_config.ppi_config.id = file

            launch_multiprocess_propgation(singular_config)

            for i in range(len(singular_config.propagations)):
                singular_config.propagations[i].id = current_id_to_original_id[i]










