import ndex2.client


from networks.loaders.base import BaseNetworkLoader
from generic_utils.constants import SpeciesIDs, NodeAttrs
from generic_utils.data_handlers.extractors import HSapiensExtractor
from generic_utils.data_handlers.translators import GeneinfoToEntrezID
from propagation.classes import PropagationNetwork, PropagationContainer


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


class HSapiensNetworkLoader(PPINetworkLoader):
    HUMAN_SPECIES_ID = SpeciesIDs.HUMAN.value

    def __init__(self, network_soruce_path):
        self.network_source_path = network_soruce_path
        self.network_extractor = HSapiensExtractor(self.network_source_path)

    def load(self, *args, **kwargs):
        network = PropagationNetwork()
        edge_triplets = self.network_extractor.extract()
        for edge_triplet in edge_triplets:
            source_node, target_node, edge_weight = self.network_extractor.unpack_triplet(edge_triplet)
            for node_id in [source_node, target_node]:
                if node_id not in network.nodes:
                    network.add_node(node_id, **self._new_node_attrs(self.HUMAN_SPECIES_ID))
            network.add_edge(source_node, target_node, **self._new_edge_attrs(edge_weight))
        return network

    @staticmethod
    def record_network(network, file_path):
        edges = sorted(network.edges(data=True))
        with open(file_path, 'w') as handler:
            for edge_data in edges:
                source = min(edge_data[0], edge_data[1])
                target = max(edge_data[0], edge_data[1])
                weight = edge_data[2][PropagationNetwork.EDGE_WEIGHT]
                handler.write(f"{source}\t{target}\t{weight}\n")


class CombinedHumanCovidNetworkLoader(HSapiensNetworkLoader):
    CORONAVIRUS_SPECIES_ID = SpeciesIDs.CORONAVIRUS.value
    NDEX_WEIGHT_KEY = "MIST"
    MERGED_COVID_NODE_NAME = "covid"

    def __init__(self, human_network_path: str, covid_human_ppi_path: str, translation_file: str,
                 merge_covid=True, covid_to_human_edge_weights=None):
        super().__init__(human_network_path)
        self.covid_human_ppi_path = covid_human_ppi_path
        self.translator = GeneinfoToEntrezID(translation_file)
        self.merge_covid = merge_covid
        self.covid_to_human_edge_weights = covid_to_human_edge_weights

    def load(self, *args, **kwargs):
        # initialize network as human only, then add in coronavirus
        network = super().load()
        covid_to_human_network = self._load_ndex()
        self._merge_covid_to_human(network, covid_to_human_network)
        return network

    def _merge_covid_to_human(self, human_network, covid_to_human_network):
        for edge in covid_to_human_network.edges(data=True):
            data = edge[2]
            if self.NDEX_WEIGHT_KEY not in data:
                continue
            source_symbol, target_symbol = data["name"].lower().split(' (interacts with) ')
            if self.merge_covid:
                source_symbol = self.MERGED_COVID_NODE_NAME
            target = str(self.translator.translate(target_symbol.upper()))
            if source_symbol not in human_network.nodes:
                human_network.add_node(source_symbol, **self._new_node_attrs(SpeciesIDs.CORONAVIRUS.value))
            # if target not in human_network.nodes:
            #     human_network.add_node(target, **self._new_node_attrs(SpeciesIDs.HUMAN.value))
            edge_weight = self.covid_to_human_edge_weights or float(data[self.NDEX_WEIGHT_KEY])
            human_network.add_edge(source_symbol, target, **self._new_edge_attrs(edge_weight))

    def _load_ndex(self):
        raw_network = ndex2.create_nice_cx_from_file(self.covid_human_ppi_path)
        return raw_network.to_networkx(mode="default")
