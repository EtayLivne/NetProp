from abc import ABCMeta, abstractmethod


class BaseNetworkLoader(metaclass=ABCMeta):
    def __init__(self, network_source_path):
        self.network_source_path = network_source_path

    @abstractmethod
    def load_network(self, **kwargs):
        """
        :param network_id: an id to be given to the generated network
        :return: an instance of PropagationGraph
        """
        raise NotImplementedError
