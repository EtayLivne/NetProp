from abc import ABCMeta, abstractmethod


class AbstractDataManager(metaclass=ABCMeta):
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def get_data(self, raw=False):
        pass