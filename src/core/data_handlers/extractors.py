import typing
from abc import ABCMeta, abstractmethod
from os import path as os_path


class AbstractDataExtractor(metaclass=ABCMeta):
    @abstractmethod
    def extract(self):
        raise NotImplemented


class BaseFileDataExtractor(AbstractDataExtractor, metaclass=ABCMeta):
    def __init__(self, file_path=None):
        self.file_path = file_path

    def _validate_path(self):
        if not os_path.isfile(self.file_path):
            raise FileNotFoundError

    @abstractmethod
    def _extract(self):
        raise NotImplemented

    def extract(self, file_path=None):
        if file_path:
            self.file_path = file_path
        self._validate_path()
        return self._extract()


class AbstractWebAPIDataExtractor(AbstractDataExtractor, metaclass=ABCMeta):
    pass
