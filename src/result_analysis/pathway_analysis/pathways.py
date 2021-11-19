import scipy.stats as sp_stats
from math import sqrt
import os
import json
from common.network_loaders import HSapeinsNetworkLoader
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from pydantic import BaseModel
from typing import Set, List, Optional
from core.data_handlers.extractors import BaseFileDataExtractor


class Pathway(BaseModel):
    tags: Optional[Set[str]]
    genes: Set[str]


class PathwayExtractor(BaseFileDataExtractor):
    def _extract(self):
        pathway_manager = PathwayManager
        with open(self.file_path, 'r') as handler:
            for line in handler.readlines():
                as_list = line.split()
                pathway_manager.add(as_list[0], set(as_list[2:]))

        return pathway_manager


class PathwayManager:
    def __init__(self):
        self._pathways = dict()

    @staticmethod
    def from_file_path(file_path):
        return PathwayExtractor(file_path).extract()

    def add(self, pathway_name: str, pathway_genes: List[str]):
        self._pathways[pathway_name] = Pathway(genes=pathway_genes)

    def get_genes(self, pathway_name, default=None):
        res = self._get(pathway_name, default=default)
        if isinstance(res, Pathway):
            return res.genes
        return res

    def get_tags(self, pathway_name, default=None):
        res = self._get(pathway_name, default=default)
        if isinstance(res, Pathway):
            return res.tags
        return res

    def tag(self, pathway_name, tag):
        try:
            pathway = self._pathways[pathway_name]
        except KeyError:
            raise ValueError(f"no pathway named {pathway_name}")
        pathway.tags.add(tag)

    def pathways(self, blacklist=None, whitelist=None):
        blacklist = blacklist or set()
        whitelist = whitelist or set()
        for pathway_name in self._pathways:
            tags = pathway_name.tags or set()
            if whitelist.issubset(tags) and not blacklist.intersection(tags):
                yield pathway_name, self._pathways[pathway_name]

    def _get(self, pathway_name, default=None):
        return self._pathways.get(pathway_name, default)


def tag_unknown_genes(pathways_manager, nodes):
    for pathway_name, pathway in pathways_manager.pathways():
        if [gene for gene in pathway.genes if gene not in nodes]:
            pathways_manager.tag(pathway_name, "unknown_genes")