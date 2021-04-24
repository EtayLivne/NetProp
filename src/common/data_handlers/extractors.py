import csv
import json
import yaml
import ndex2.client
from mygene import MyGeneInfo
from core.data_handlers.extractors import BaseFileDataExtractor, AbstractWebAPIDataExtractor


class HSapiensExtractor(BaseFileDataExtractor):
    def _extract(self):
        with open(self.file_path, 'r') as handler:
            return [line.split()[:3] for line in handler.readlines()]


class CSVExtractor(BaseFileDataExtractor):
    def _extract(self):
        with open(self.file_path, 'r') as handler:
            return {i: row_dict for i, row_dict in enumerate(csv.DictReader(handler))}


class JsonExtractor(BaseFileDataExtractor):
    def _extract(self):
        with open(self.file_path, 'r') as handler:
            return json.load(handler)

    @staticmethod
    def dump(data, file_path):
        with open(file_path, 'w') as handler:
            json.dump(data, handler, indent=4)


class YamlExtractor(BaseFileDataExtractor):
    def _extract(self):
        with open(self.file_path, 'r') as handler:
            return yaml.load(handler)

    @staticmethod
    def dump(data, file_path):
        with open(file_path, 'w') as handler:
            yaml.dump(data, handler, indent=4)


class NDEXExtractor(AbstractWebAPIDataExtractor):
    def __init__(self, server_url=None, network_id=None):
        self._server_url = server_url
        self._network_id = network_id

    def extract(self, server_url=None, network_id=None, as_nx=True):
        self._server_url = server_url or self._server_url
        self._network_id = network_id or self._network_id
        if not (self._server_url or self._network_id):
            raise ValueError("server and network ID must be specified in order to extract NDEX data")
        try:
            nice_cx_from_server = ndex2.create_nice_cx_from_server(server=self._server_url,
                                                                   uuid=self._network_id)
        except:
            raise ValueError(f'attempt to query network {self._network_id} from NDEX server {self._server_url} failed')

        return nice_cx_from_server.to_networkx(mode="default") if as_nx else nice_cx_from_server


class GeneInfoExtractor(AbstractWebAPIDataExtractor):
    def __init__(self):
        pass

    def extract(self, gene_names_list):
        try:
            ncbi_query = MyGeneInfo().querymany(gene_names_list, scopes="symbol", fields=["entrezgene", "symbol"],
                                                species="human")
        except:
            raise Exception("Failed to query MyGeneInfo")
        return [result for result in ncbi_query]

