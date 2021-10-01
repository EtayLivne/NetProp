from pydantic import BaseModel, validator
from typing import Optional, Set, List
from enum import Enum

class ProteinIDs(str, Enum):
    ENTREZGENE = "entrezgene"


class HaltConditionOptions(str, Enum):
    NUMBER_OF_ROUNDS = "number_of_rounds"
    GAP_THRESHOLD = "gap_threshold"


class PPIModel(BaseModel):
    network_id: str = "ppi"
    source_file: str
    network_loader: str
    network_loader_input_dict: Optional[dict] = None
    directed = False
    protein_id_class: ProteinIDs


class PropagationSourceModel(BaseModel):
    id: str
    source_of: Set[str]


class PropagationTargetModel(BaseModel):
    id: str
    target_of: Set[str]


class HaltConditionOptionModel(BaseModel):
    type: HaltConditionOptions
    number_of_rounds: Optional[int]
    gap_threshold: Optional[float]

    @validator("number_of_rounds", "threshold_gap")
    def is_number_of_rounds_condition(cls, v, values, **kwargs):
        if v is not None and "type" in values and values["type"] != v:
            raise ValueError("in halt condition of type {} field {} is not applicable".
                             format(values["type"], v))


class PropagationParametersModel(BaseModel):
    prior_set: Set[PropagationSourceModel]
    target_set: Set[PropagationTargetModel]
    prior_set_confidence: float
    halt_conditions: List[HaltConditionOptionModel]
    propagation_id: Optional[str] = None

    @validator('prior_set_confidence')
    def confidence_is_valid_probability(cls, v):
        if v < 0 or v > 1:
            raise ValueError("the confidence in a prior set must be a probability between 0 and 1")
        return v


class ConfigModel(BaseModel):
    ppi_config: PPIModel
    propagations: List[PropagationParametersModel]
    output_dir_path: str

