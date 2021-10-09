from pydantic import BaseModel, validator
from typing import Optional, Set, List, Dict
from enum import Enum
from config_models import ProteinIDs, HaltConditionOptionModel
from constants import SpeciesIDs


class PropagationNetworkModel(BaseModel):
    source_file: str
    network_id: str
    directed: bool
    suppressed_nodes: Optional[List[str]]


class PropagationNodeModel(BaseModel):
    source_of: Set[str]
    target_of: Set[str]
    liquids: Dict[str, float]


class PropagationProteinModel(PropagationNodeModel):
    id_type: ProteinIDs
    species: SpeciesIDs


class PropagationResultModel(BaseModel):
    network: PropagationNetworkModel
    nodes: List[PropagationNodeModel]
    prior_set_confidence: float
    halt_gap_conditions: List[HaltConditionOptionModel]
