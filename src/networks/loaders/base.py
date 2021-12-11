from abc import ABCMeta, abstractmethod
from models.config_models import NetworksParametersModel
from propagation.classes import PropagationNetwork
from typing import List


class BaseNetworkLoader(metaclass=ABCMeta):
    def __init__(self, *args, **kwargs):
        pass


class SingleNetworkLoader(BaseNetworkLoader, metaclass=ABCMeta):
    @abstractmethod
    def load(self, *args, **kwargs) -> PropagationNetwork:
        """
        :param network_id: an id to be given to the generated network
        :return: an instance of PropagationGraph
        """
        raise NotImplementedError


class MultiNetworkLoader(BaseNetworkLoader, metaclass=ABCMeta):
    @abstractmethod
    def load(self, *args, **kwargs) -> List[NetworksParametersModel]:
        """
        :param network_id: an id to be given to the generated network
        :return: an instance of PropagationGraph
        """
        raise NotImplementedError