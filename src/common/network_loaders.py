from network_loader import BaseNetworkLoader
from common.data_handlers.extractors import HSapiensExtractor
from common.data_handlers.translators import GeneinfoToEntrezID
from propagater import PropagationNetwork, PropagationContainer
from constants import SpeciesIDs, NodeAttrs
import ndex2.client
import networkx as nx

class PPINetworkLoader(BaseNetworkLoader):
    SPECIES_ID_KEY = NodeAttrs.SPECIES_ID.value

    @classmethod
    def _new_node_attrs(cls, species):
        return {
            PropagationNetwork.CONTAINER_KEY: PropagationContainer(),
            cls.SPECIES_ID_KEY: species
        }

    @classmethod
    def _new_edge_attrs(cls, weight):
        return {
            PropagationNetwork.EDGE_WEIGHT: weight
        }

class HSapeinsNetworkLoader(PPINetworkLoader):
    HUMAN_SPECIES_ID = SpeciesIDs.HUMAN.value

    def __init__(self, network_soruce_path):
        super().__init__(network_soruce_path)
        self.network_extractor = HSapiensExtractor(self.network_source_path)

    def load_network(self, **kwargs):
        network = PropagationNetwork()
        edge_triplets = self.network_extractor.extract()
        for edge_triplet in edge_triplets:
            source_node, target_node, edge_weight = self.network_extractor.unpack_triplet(edge_triplet)
            for node_id in [source_node, target_node]:
                if node_id not in network.nodes:
                    network.add_node(node_id, **self._new_node_attrs(self.HUMAN_SPECIES_ID))
            network.add_edge(source_node, target_node, **self._new_edge_attrs(edge_weight))
        return network


class CovidToHumanNetworkLoader(PPINetworkLoader):
    _DEFAULT_PATH_TO_ID_TRANSLATION_FILE = r"C:\studies\code\NetProp\data\gene_symbol_to_entrezgeneid_new.json"

    def load_network(self, **kwargs):
        raw_network = ndex2.create_nice_cx_from_file(self.network_source_path)
        nx_network = raw_network.to_networkx(mode="default")
        if kwargs.get("normalize", False):
            normalized_network = nx.Graph()
            path_to_translation_file = kwargs.get("translation_file_path", self._DEFAULT_PATH_TO_ID_TRANSLATION_FILE)
            symbol_to_id_translator = GeneinfoToEntrezID(path_to_translation_file)
            debug_count = 0
            for edge in nx_network.edges(data=True):
                data = edge[2]
                if "MIST" not in data:
                    continue
                source_symbol, target_symbol = data["name"].lower().split(' (interacts with) ')
                try:
                    target = str(symbol_to_id_translator.translate(target_symbol))
                except ValueError:
                    print(f"could not add edge for target symbol {target_symbol}")
                if source_symbol not in normalized_network.nodes:
                    debug_count += 1
                    normalized_network.add_node(source_symbol, **self._new_node_attrs(SpeciesIDs.CORONAVIRUS.value))
                if target not in normalized_network.nodes:
                    normalized_network.add_node(target, **self._new_node_attrs(SpeciesIDs.HUMAN.value))
                normalized_network.add_edge(source_symbol, target, **self._new_edge_attrs(float(data['MIST'])))
            return normalized_network
        else:
            return nx_network


class HumanCovidHybridNetworkLoader(HSapeinsNetworkLoader):
    CORONAVIRUS_SPECIES_ID = SpeciesIDs.CORONAVIRUS.value
    _DEFAULT_PATH_TO_COVID_PPI_FILE = r"C:\studies\code\NetProp\data\ndex_covid_human_ppi.cx"

    def load_network(self, **kwargs):
        # initialize network as human only, then add in coronavirus
        hybrid_network = super().load_network()
        covid_ppi_path = kwargs.get("covid_ppi_path", self._DEFAULT_PATH_TO_COVID_PPI_FILE)
        covid_to_human_network = CovidToHumanNetworkLoader(covid_ppi_path).load_network(normalize=True)
        hybrid_network.add_nodes_from([(node, data) for node, data in covid_to_human_network.nodes(data=True) if
                                       data[NodeAttrs.SPECIES_ID.value] == SpeciesIDs.CORONAVIRUS.value])
        hybrid_network.add_edges_from(covid_to_human_network.edges(data=True))
        return hybrid_network
