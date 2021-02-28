import networkx as nx

from core.data_handlers.managers import AbstractDataManager
from common.data_handlers.extractors import HSapiensExtractor, CSVExtractor, JsonExtractor,  NDEXExtractor,\
                                            GeneInfoExtractor


class HSapiensManager(AbstractDataManager):
    def __init__(self, file_path):
        self.extractor = HSapiensExtractor(file_path=file_path)

    def _get_raw_data(self):
        return self.extractor.extract()

    def get_data(self, raw=False):
        raw_data = self._get_raw_data()
        if raw:
            return raw_data

        source_node_index = 0
        target_node_index = 1
        edge_weight_index = 2

        # TODO figure out graph structure than construct accordingly
        graph = nx.DiGraph()
        discovered_nodes_map = dict()
        for triplet in raw_data:
            for node_index in (triplet[source_node_index], triplet[target_node_index]):
                if node_index not in discovered_nodes_map:
                    discovered_nodes_map[node_index] = Protein(id=int(node_index))

            graph.add_edge(discovered_nodes_map[triplet[source_node_index]],
                           discovered_nodes_map[triplet[target_node_index]],
                           weight=float(triplet[edge_weight_index]))
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

