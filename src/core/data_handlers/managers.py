from abc import ABCMeta, abstractmethod


class AbstractDataManager(metaclass=ABCMeta):
    @abstractmethod
    def get_data(self, raw=False):
        raise NotImplemented
