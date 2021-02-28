from abc import ABCMeta, abstractmethod
from math import floor

class AbstractDataTranslator(metaclass=ABCMeta):
    @abstractmethod
    def translate(self, data):
        raise NotImplemented

    @abstractmethod
    def _init_translation_map(self):
        raise NotImplemented


class AbstractTranslationMap(metaclass=ABCMeta):
    @abstractmethod
    def translate(self, key):
        raise NotImplemented

    @abstractmethod
    def reverse_translate(self, key):
        raise NotImplemented

    def translate_many(self, keys):
        return {key: self.translate(key) for key in keys}


class BaseTranslationMap(metaclass=ABCMeta):
    def __init__(self, translation_dict):
        self._map = dict()
        self._reverse_map = dict()
        self._reverse_cache_max_size_ratio = 0.1

    def _insert_to_reverse_map(self, to_insert=dict(), force=False):
        if not isinstance(to_insert, dict):
            raise ValueError("only dicts can be merged into a reverse translation map")
        if not force and {k for k in to_insert if k in self._reverse_map}:
            raise ValueError("attempt to override elements in reverse translation map!")

        self._reverse_map.update(to_insert)
        excess_items = len(self._reverse_map) - floor(len(self._map) * self._reverse_cache_max_size_ratio)
        keys_to_remove = [k for k in self._reverse_map if k not in to_insert][:max(0, excess_items)]
        for k in keys_to_remove:
            self._reverse_map.pop(k)

    def resize_reverse_cache(self, new_size_ratio):
        self._reverse_cache_max_size_ratio = new_size_ratio
        while len(self._reverse_map) > self._reverse_cache_max_size_ratio:
            self.insert_to_reverse_map(None)

    def translate(self, key):
        try:
            return self._map[key]
        except KeyError:
            raise ValueError(f'Could not translate {key} because it does not exist in the translation map')

    def reverse_translate(self, key):
        if key not in self._reverse_map:
            try:
                value = next(k for k in self._map if self._map[k] == key)
            except StopIteration:
                raise ValueError(f'Could not translate {key} because it does not exist in the translation map')
            self._insert_to_reverse_map({key: value})
        return self._reverse_map[key]


