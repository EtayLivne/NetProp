from typing import List
from abc import ABCMeta, abstractmethod

from ..propagation_network import PropagationNetwork
from netprop.models.config_models import NetworksParametersModel

loader_registry = dict()

class BaseNetworkLoader(metaclass=ABCMeta):

    def __init__(self, *args, **kwargs):
        pass


    def __init_subclass__(cls, **kwargs):
        loader_registry[cls.__name__] = cls


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