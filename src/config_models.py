from pydantic import BaseModel, validator
from typing import Optional, Set, List, ClassVar
from enum import Enum

class ProteinIDs(str, Enum):
    ENTREZGENE = "entrezgene"


class HaltConditionOptions(str, Enum):
    MAX_STEPS = "max_steps"
    MIN_GAP = "min_gap"


class PPIModel(BaseModel):
    source_file_key: ClassVar[str] = "source_file"
    directed_key: ClassVar[str] = "directed"
    protein_id_class_key: ClassVar[str] = "protein_id_class"
    network_id_key: ClassVar[str] = "network_id"

    network_id: str = "ppi"
    source_file: str
    network_loader: str
    network_loader_input_dict: Optional[dict] = None
    directed = False
    protein_id_class: ProteinIDs
    network_id: Optional[str]


class PropagationSourceModel(BaseModel):
    id: str
    source_of: Set[str]


class PropagationTargetModel(BaseModel):
    id: str
    target_of: Set[str]


class HaltConditionOptionModel(BaseModel):
    condition_type: HaltConditionOptions
    max_steps: Optional[int]
    min_gap: Optional[float]

    @validator("max_steps", "min_gap", allow_reuse=True)
    def condition_matches_type(cls, v, values, **kwargs):
        if v is not None and "condition_type" in values and values["condition_type"] != v:
            raise ValueError("in halt condition of type {} field {} is not applicable".format(values["type"], v))


class PropagationParametersModel(BaseModel):
    prior_set: Set[PropagationSourceModel]
    target_set: Set[PropagationTargetModel]
    suppressed_set: Set[str]
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

