from pydantic import BaseModel
from typing import Optional, Set, List, Dict
from config_models import ProteinIDs, HaltConditionOptionModel
from constants import SpeciesIDs


class PropagationNetworkModel(BaseModel):
    source_file: str
    network_id: str
    directed: bool
    suppressed_nodes: Optional[List[str]]


class PropagationNodeModel(BaseModel):
    id: str
    source_of: Set[str]
    target_of: Set[str]
    liquids: Dict[str, float]

    class Config:
        json_encoders = {set: list}

    def __hash__(self):
        return int(self.id)

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()


class PropagationProteinModel(PropagationNodeModel):
    id_type: ProteinIDs
    species: SpeciesIDs

    class Config(PropagationNodeModel.Config):
        use_enum_values = True


class PropagationResultModel(BaseModel):
    id: str
    network: PropagationNetworkModel
    nodes: List[PropagationNodeModel]
    prior_set_confidence: float
    halt_conditions: HaltConditionOptionModel

