from core.data_handlers.translators import AbstractDataTranslator, BaseTranslationMap
from common.data_handlers.extractors import JsonExtractor, CSVExtractor


class GeneinfoToEntrezID(AbstractDataTranslator):
    def __init__(self, path_to_translation_file):
        self._extractor = JsonExtractor(path_to_translation_file)
        self._translation_map = self._init_translation_map()

    def _init_translation_map(self):
        return BaseTranslationMap(self._extractor.extract())

    def translate(self, key):
        return self._translation_map.translate(key)

    def reverse_translate(self, key):
        return self._translation_map.reverse_translate(key)


class ProteinRoles(AbstractDataTranslator):
    def __init__(self, path_to_translation_file):
        self.extractor = CSVExtractor(path_to_translation_file)
        self.translation_map = self._init_translation_map()

    def _init_translation_map(self):
        translation_dict = {row['protein']: {k for k in row if row[k]} for row in self.extractor.extract()}
        return BaseTranslationMap(translation_dict)

    def translate(self, key):
        return self.translation_map.translate(key)


class DrugbankIDToDrugNames(AbstractDataTranslator):
    DRUGBANK_ID_COL = "DrugBank ID"
    ACCESSION_NUMBER_COL = "Accession Numbers"
    COMMON_NAME_COL = "Common name"
    CAS_COL = "CAS"
    UNII_COL = "UNII"
    SYNONYMS_COL = "Synonyms"
    INCHI_COL = "Standard InChI Key"

    def __init__(self, path_to_translation_file):
        self.extractor = CSVExtractor(path_to_translation_file)
        self.translation_map = self._init_translation_map()

    def _init_translation_map(self):
        raw_data = self.extractor.extract()
        translators = {k: dict() for k in raw_data[0] if k != self.DRUGBANK_ID_COL}
        for row in raw_data:
            try:
                drugbank_id = row.pop(self.DRUGBANK_ID_COL)
            except KeyError:
                raise Exception(f'drugbank id to drug names translation map source file format error '
                                f'(no column named {self.DRUGBANK_ID_COL})')
            for k, v in row.items():
                translators[k].update({drugbank_id: v})

        self.translation_map = {k: BaseTranslationMap(translators[k]) for k in translators}

    def translate(self, data, accession_numbers=False, common_name=False, cas=False, unii=False, synonyms=False, inchi_key=False):
        return_set = set()
        no_filter = not (accession_numbers or common_name or cas or unii or inchi_key)

        if no_filter or accession_numbers:
            return_set.add(self.translation_map[self.ACCESSION_NUMBER_COL].translate())
        if no_filter or common_name:
            return_set.add(self.translation_map[self.COMMON_NAME_COL].translate())
        if no_filter or cas:
            return_set.add(self.translation_map[self.CAS_COL].translate())
        if no_filter or unii:
            return_set.add(self.translation_map[self.UNII_COL].translate())
        if no_filter or synonyms:
            return_set.add(self.translation_map[self.SYNONYMS_COL].translate())
        if no_filter or inchi_key:
            return_set.add(self.translation_map[self.INCHI_COL].translate())

        if len(return_set) != 1:
            return return_set
        return return_set.pop()

    # TODO reverse translate?
