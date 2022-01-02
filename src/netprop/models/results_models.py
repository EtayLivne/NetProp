from pydantic import BaseModel
from typing import Optional, Set, List, Dict

from .config_models import HaltConditionOptionModel


class PropagationNetworkModel(BaseModel):
    source_file: str
    network_id: str
    directed: bool
    suppressed_nodes: Optional[List[str]]


class PropagationNodeModel(BaseModel):
    id: str
    source_of: Set[str]
    # target_of: Set[str] #TODO delete if /when we retire target sets
    rank: Optional[int]
    liquids: Dict[str, float]
    metadata: Dict

    class Config:
        json_encoders = {set: list}
        use_enum_values = True

    def __hash__(self):
        return int(self.id)

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()


class PropagationResultModel(BaseModel):
    id: str
    network: dict
    nodes: Dict[str, PropagationNodeModel]
    prior_set_confidence: float
    halt_conditions: HaltConditionOptionModel

