from pydantic import BaseModel, validator, root_validator
from typing import Any
from typing import List, Dict


class NetpropNodeModel(BaseModel):
    id: str
    data: Dict[str, Any]


class NetpropEdgeModel(BaseModel):
    source: str
    target: str
    weight: float

class NetpropNetworkModel(BaseModel):
    nodes: List[NetpropNodeModel]
    edges: List[NetpropEdgeModel]
    data: Dict[str, Any]

    @root_validator
    def edges_only_between_known_nodes(cls, values):
        nodes, edges = values.get('nodes'), values.get('edges')
        node_ids = {node.id for node in nodes}
        edge_ids = {e.source for e in edges} | {e.target for e in edges}
        bad = edge_ids - node_ids
        if bad:
            raise ValueError(f"The following nodes are referenced in edges but do not exist in the graph: {bad}")
        return values