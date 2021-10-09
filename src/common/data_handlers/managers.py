import networkx as nx

from common.data_classes.data_classes import Protein
from core.data_handlers.managers import AbstractDataManager
from common.data_handlers.extractors import HSapiensExtractor, CSVExtractor, NDEXExtractor, GeneInfoExtractor,\
                                            JsonExtractor
from common.data_classes.data_classes import HUMAN_SPECIES_NAME
from propagater import PropagationNetwork, PropagationResult


class HSapiensManager(AbstractDataManager):
    def __init__(self, file_path):
        self.extractor = HSapiensExtractor(file_path=file_path)

    def _get_raw_data(self):
        return self.extractor.extract()

    def get_data(self, raw=False):
        raw_data = self._get_raw_data()
        if raw:
            return raw_data

        human_graph = PropagationNetwork()
        for triplet in raw_data:
            source_node, target_node, edge_weight = int(triplet[0]), int(triplet[1]), float(triplet[2])
            human_graph.add_edge(source_node, target_node, weight=edge_weight)
            for node_id in [source_node, target_node]:
                if not human_graph.nodes[node_id]:
                    human_graph.nodes[node_id][human_graph.CONTAINER_KEY] = \
                        Protein(id=node_id, species=HUMAN_SPECIES_NAME)

        return human_graph


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


class PropagationResultsManager(AbstractDataManager):
    def __init__(self, file_path=None, propagation_results=None):
        self.extractor = JsonExtractor(file_path=file_path)
        self.propagation_results = propagation_results or []

    def reload_from_file(self, file_path=None):
        self.propagation_results = self.extractor.extract(file_path=file_path)

    def _get_raw_data(self, reload_from_file=False):
        if reload_from_file:
            self.reload_from_file()
        return self.propagation_results

    def get_data(self, raw=False):
        return self.propagation_results

    def dump_to_file(self, file_path):
        self.extractor.dump(self.propagation_results, file_path)