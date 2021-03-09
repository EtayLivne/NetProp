import networkx as nx

from common.data_classes.data_classes import Protein
from core.data_handlers.managers import AbstractDataManager
from common.data_handlers.extractors import HSapiensExtractor, CSVExtractor, NDEXExtractor, GeneInfoExtractor



class HSapiensManager(AbstractDataManager):
    def __init__(self, file_path):
        self.extractor = HSapiensExtractor(file_path=file_path)

    def _get_raw_data(self):
        return self.extractor.extract()

    def get_data(self, raw=False):
        raw_data = self._get_raw_data()
        if raw:
            return raw_data

        graph = nx.DiGraph()
        for triplet in raw_data:
            source_node, target_node, edge_weight = triplet
            graph.add_edge(source_node, target_node, weight=float(edge_weight))
            for node_id in [source_node, target_node]:
                if not g.nodes[node_id]:
                    g.nodes[node_id]["protein_data"] = Protein(id=int(node_id))

        return graph


class NDEXManager(AbstractDataManager):
    def __init__(self, server_url):
        self.server_url = server_url
        self.extractor = NDEXExtractor(server_url=server_url)

    def get_data(self, network_id, raw=False):
        return self.extractor.extract(network_id=network_id, as_nx=raw)


class GeneInfoManager(AbstractDataManager):
    def __init__(self):
        self.extractor = GeneInfoExtractor()

    def get_data(self, gene_names):
        return self.extractor.extract(list(gene_names))


class DrugbankTargetsManager(AbstractDataManager):
    def __init__(self, file_path):
        self.extractor = CSVExtractor(file_path=file_path)

    def _get_raw_data(self):
        return self.extractor.extract()

    def get_data(self, raw=False):
        raw_data = self._get_raw_data()
        if raw:
            return raw_data
        # TODO finish after defining data class for drugbank target information

