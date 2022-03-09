import json
import ndex2.client
import networkx as nx


from .base import SingleNetworkLoader
from netprop.models import NetpropNetworkModel, NetpropNodeModel, NetpropEdgeModel
from netprop.generic_utils.constants import SpeciesIDs, NodeAttrs
from netprop.generic_utils.data_handlers.extractors import HSapiensExtractor
from netprop.generic_utils.data_handlers.translators import GeneinfoToEntrezID
from netprop.propagation.classes import PropagationNetwork, PropagationContainer


class PPINetworkLoader(SingleNetworkLoader):
    SPECIES_ID_KEY = NodeAttrs.SPECIES_ID.value

    @classmethod
    def _new_node_attrs(cls, species):
        return {
           "source_of": [],
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
        network = nx.Graph()
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


class MetaCovidHumanLoader(HSapiensNetworkLoader):
    MERGED_COVID_NODE_NAME = "covid"
    CONFIDENCE_SCORES = {
        2: 0.8,
        3: 0.85,
        4: 0.9,
        5: 0.95,
        6: 1
    }

    def __init__(self, human_network_path: str, covid_to_human_path: str, merged_covid=True):
        super().__init__(human_network_path)
        self.covid_to_human_path = covid_to_human_path
        self.merged_covid = merged_covid

    def load(self, *args, **kwargs):
        network = super().load()
        with open(self.covid_to_human_path, 'r') as handler:
            cov_to_human_ppi = json.load(handler)
        for cov_protein, interacting_human_proteins in cov_to_human_ppi.items():
            source = self.MERGED_COVID_NODE_NAME if self.merged_covid else cov_protein
            if source not in network:
                network.add_node(source, **self._new_node_attrs(SpeciesIDs.CORONAVIRUS.value))
            for human_protein, num_paper_appearances in interacting_human_proteins.items():
                if num_paper_appearances < 2:
                    continue
                target = human_protein
                if target not in network:
                    print(f"cannot add covid interaction with {human_protein} - it isn't in the network")
                    continue
                confidence = self._calc_confidence(num_paper_appearances)
                network.add_edge(source, target, weight=confidence)
        return network

    def _calc_confidence(self, num_paper_appearances):
        return self.CONFIDENCE_SCORES[num_paper_appearances]



class NetpropNetwork(PPINetworkLoader):
    def __init__(self, network_path: str):
        self.network_path = network_path

    def load(self, *args, **kwargs):
        network_data = NetpropNetworkModel.parse_file(self.network_path)
        network = PropagationNetwork()

        network.add_edges_from([(n.id, n.data) for n in network_data.nodes])
        network.add_weighted_edges_from([e.source, e.target, e.weight] for e in network_data.edges)
        for k, v in network_data.data.items():
            network.graph[k] = v

    @staticmethod
    def record_network(network: nx.Graph, file_path: str):
        nodes = [NetpropNodeModel(id=n, data=data) for n, data in network.nodes(data=True)]
        edges = [NetpropEdgeModel(source=e[0], target=e[1], weight=e[2]) for e in network.edges.data("weight")]
        data = network.graph
        model = NetpropNetworkModel(nodes=nodes, edges=edges, data=data)
        with open(file_path, 'w') as handler:
            json.dump(model.dict(), handler, indent=4)