from pydantic import BaseModel
from typing import Optional, Set, List, Dict, Union
from pandas import DataFrame
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

    def prop_scroes_as_series(self, by_liquids: Union[str, list[str]]="info", sort: bool=False):
        if not isinstance(by_liquids, list):
            by_liquids = [by_liquids]

        n_list = list(self.nodes.keys())
        df_columns = [n_list] + [[self.nodes[n].liquids.get(liquid, 0) for n in n_list] for liquid in by_liquids]
        df_column_names = ["nodes"] + [self.id + liquid for liquid in by_liquids]

        dict_form = {df_column_names[i]: df_columns[i] for i in range(len(df_columns))}
        df = DataFrame(dict_form)
        if sort:
            df.sort_values(ascending=False, inplace=True)
        return df
