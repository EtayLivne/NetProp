from pydantic import BaseModel, validator, Field
from typing import Optional, Set, List, ClassVar
from propagater import Propagater
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
    id_key: ClassVar[str] = "network_id"

    source_file: str
    loader_class: str
    loader_input_dict: Optional[dict] = None
    directed = False
    protein_id_class: ProteinIDs
    id: Optional[str]


class BasePropagationNodeModel(BaseModel):
    id: str
    def __hash__(self):
        return int(self.id)

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    class Config:
        json_encoders = {set: list}


class PropagationSourceModel(BasePropagationNodeModel):
    source_of_key: ClassVar[str] = "source_of"

    source_of: Set[str]


class PropagationTargetModel(BasePropagationNodeModel):
    target_of_key: ClassVar[str] = "target_of"

    target_of: Set[str]


class HaltConditionOptionModel(BaseModel):
    #conditio_type: HaltConditionOptions
    max_steps: Optional[int]
    min_gap: Optional[float]

    # @validator("max_steps", "min_gap", allow_reuse=True)
    # def condition_matches_type(cls, v, values, **kwargs):
    #     if v is not None and "condition_type" in values and values["condition_type"] != v:
    #         raise ValueError("in halt condition of type {} field {} is not applicable".format(values["type"], v))


class PropagationParametersModel(BaseModel):
    prior_set: Optional[List[PropagationSourceModel]]
    target_set: List[PropagationTargetModel] = Field(default_factory=set)
    suppressed_set: Set[str] = Field(default_factory=set)
    prior_set_confidence: Optional[float]
    halt_conditions: Optional[HaltConditionOptionModel]
    method: str = Propagater.ITERATIVE
    id: Optional[str] = None
    output_dir_path: Optional[str]

    # this model has all optional fields intentionally so that they may be filled after creation.
    # this method validates that all fields are ready
    def validate_completeness(self):
        return self.prior_set and self.prior_set_confidence and self.halt_conditions and self.output_dir_path

    @validator('prior_set_confidence')
    def confidence_is_valid_probability(cls, v):
        if v and (v < 0 or v > 1):
            raise ValueError("the confidence in a prior set must be a probability between 0 and 1")
        return v


class ConfigModel(BaseModel):
    ppi_config: PPIModel
    global_propagation_params: Optional[PropagationParametersModel]
    propagations: List[PropagationParametersModel]


