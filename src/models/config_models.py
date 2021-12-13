from pydantic import BaseModel, validator, Field
from typing import Optional, Set, List, ClassVar, Union
from propagation.classes import Propagater
from enum import Enum

class ProteinIDs(str, Enum):
    ENTREZGENE = "entrezgene"


class PropagationSourceModel(BaseModel):
    source_of_key: ClassVar[str] = "source_of"

    id: str
    source_of: List[str]

    def __hash__(self):
        return int(self.id)

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    class Config:
        json_encoders = {set: list}


class SuppressedSetModel(BaseModel):
    id: Optional[str]
    nodes: List[str] = Field(default_factory=list)


class PriorSetModel(BaseModel):
    id: Optional[str]
    nodes: List[PropagationSourceModel]
    confidence: Optional[float]

    @validator('confidence')
    def confidence_is_valid_probability(cls, v):
        if v and (v < 0 or v > 1):
            raise ValueError("the confidence in a prior set must be a probability between 0 and 1")
        return v


class HaltConditionOptionModel(BaseModel):
    max_steps: Optional[int]
    min_gap: Optional[float]


class NetworkInclusionParametersModel(BaseModel):
    path: Union[str, List[str]]
    exclude_files: Optional[List[str]]
    include_files: Optional[List[str]]
    exclude_ids: Optional[List[str]]
    include_ids: Optional[List[str]]


class NetworksParametersModel(BaseModel):
    id_key: ClassVar[str] = "network_id"

    id: Optional[str]
    multi: bool = False
    loader_class: str
    init_args: List[str] = Field(default_factory=list)
    init_kwargs: dict = Field(default_factory=dict)
    metadata: dict = Field(default_factory=dict)


class ConfigModel(BaseModel):
    networks: NetworkInclusionParametersModel
    prior_set: Union[PriorSetModel, List[PriorSetModel]]
    suppressed_set: Optional[Union[SuppressedSetModel, List[SuppressedSetModel]]]
    halt_conditions: HaltConditionOptionModel
    method: str = Propagater.ITERATIVE
    id: Optional[str] = None
    output_dir_path: Optional[str]

