from network_loader import BaseNetworkLoader
from common.data_handlers.extractors import HSapiensExtractor
from propagater import PropagationNetwork, PropagationContainer
from constants import SpeciesIDs, NodeAttrs


class PPINetworkLoader(BaseNetworkLoader):
    SPECIES_ID_KEY = NodeAttrs.SPECIES_ID.value


class HSapeinsNetworkLoader(PPINetworkLoader):
    HUMAN_SPECIES_ID = SpeciesIDs.HUMAN.value

    def __init__(self, network_soruce_path):
        super().__init__(network_soruce_path)
        self.network_extractor = HSapiensExtractor(self.network_source_path)

    @classmethod
    def _new_node_attrs(cls):
        return {
            PropagationNetwork.CONTAINER_KEY: PropagationContainer(),
            cls.SPECIES_ID_KEY: cls.HUMAN_SPECIES_ID
        }

    @classmethod
    def _new_edge_attrs(cls, weight):
        return {
            PropagationNetwork.EDGE_WEIGHT: weight
        }

    def load_network(self, **kwargs):
        network = PropagationNetwork()
        edge_triplets = self.network_extractor.extract()
        for edge_triplet in edge_triplets:
            source_node, target_node, edge_weight = self.network_extractor.unpack_triplet(edge_triplet)
            for node_id in [source_node, target_node]:
                if node_id not in network.nodes:
                    network.add_node(node_id, **self._new_node_attrs())
            network.add_edge(source_node, target_node, **self._new_edge_attrs(edge_weight))
        return network

