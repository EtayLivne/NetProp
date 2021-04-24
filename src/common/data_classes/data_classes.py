import os
from enum import Enum
from copy import deepcopy

from dataclasses import dataclass, field
from functools import total_ordering
from propagater import PropagationContainer

HUMAN_SPECIES_NAME = "human"
COVID_19_SPECIES_NAME = "covid"

@dataclass()
class Protein(PropagationContainer):
    species: str = HUMAN_SPECIES_NAME


@dataclass()
class DrugBankProteinTarget:
    id: str
    name: str
    gene_name: str
    GenBank_Protein_id: str
    GenBank_gene_id: str
    UniProt_id: str
    UniProt_title: str
    PDB_id: str
    GeneCard_id: str
    GenAtlas_id: str
    HGNC_id: str
    species: str
    drug_ids: list = field(default_factory=list)


@dataclass(frozen=True)
class GeneSet:
    name: str
    gene_set: set = field(default_factory=set)

    def __iter__(self):
        return iter(self.gene_set)


class RunConfig:
    class ConfKeys(Enum):
        NETWORK = "network"
        NETWORKS_DIR = "networks_dir"
        PRIOR_SET = "prior_set"
        PRIOR_SETS = "prior_sets"
        KNOCKOUT_SET = "knockout_set"
        KNOCKOUT_SETS = "knockout_sets"
        PRIOR_SET_CONFIDENCE = "prior_set_confidence"
        OUTPUT_DIR = "output_dir"

    class SetConfKeys(Enum):
        SET = "set"
        CONFIDENCE = "confidence"
        REPEATS = "repeats"
        RANDOMIZED = "randomized"

    def __init__(self, config_dict):
        if len(config_dict) != 1:
            raise ValueError("a config dict must be of the format {name: config_class}")
        self.name = config_dict.keys()[0]

        config = config_dict[self.name]
        if unknown_keys := set(config.keys()) - {e.value for e in self.ConfKeys}:
            raise ValueError(f'The following kwargs in config of run {self.name} are undefined: {unknown_keys}')
        contradictory_pairs = [(self.ConfKeys.NETWORK, self.ConfKeys.NETWORKS_DIR),
                               (self.ConfKeys.PRIOR_SET, self.ConfKeys.PRIOR_SETS),
                               (self.ConfKeys.KNOCKOUT_SET_SET, self.ConfKeys.KNOCKOUT_SETS)]
        if contradictory_pair := next(pair for pair in contradictory_pairs if pair[0] in config and pair[1] in config):
            raise ValueError(f"config of run {self.name} has contradictory kwargs pair: {contradictory_pair}")

        if network_path := config.get(self.ConfKeys.NETWORK, False):
            self.networks = [network_path]
        elif networks_dir := config.get(self.ConfKeys.NETWORKS_DIR, False):
            self.networks = [os.path.join(networks_dir, file) for file in os.listdir(networks_dir)]
        else:
            raise ValueError(f'config of run {self.name} is missing a kwarg for source network or networks dir')

        self.prior_sets = config.get(self.ConfKeys.PRIOR_SETS, None) or list(config.get(self.ConfKeys.PRIOR_SET, None))
        self.knockout_sets = config.get(self.ConfKeys.KNOCKOUT_SETS, None) or list(config.get(self.ConfKeys.KNOCKOUT_SET, None))
        if self.knockout_sets and len(self.prior_sets) != len(self.knockout_sets):
            raise ValueError(f'config of run {self.name} has prior a different number of knockout sets and prior sets')

        self.default_run_prior_set_confidence = config.get(self.ConfKeys.PRIOR_SET_CONFIDENCE, None)
        self.output_dir = config.get(self.ConfKeys.OUTPUT_DIR, os.getcwd())


class Config:
    def __init__(self, config_data):
        run_config_list = config_data if type(config_data) is list else [config_data]
        self.run_config_list = [RunConfig(deepcopy(run_config_dict)) for run_config_dict in run_config_list]





