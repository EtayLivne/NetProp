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
        if not (accession_numbers or common_name or cas or unii or synonyms or inchi_key):
            accession_numbers = common_name = cas = synonyms = inchi_key = True

        param_to_col_map = {
            accession_numbers: self.ACCESSION_NUMBER_COL,
            common_name: self.COMMON_NAME_COL,
            cas: self.CAS_COL,
            unii: self.UNII_COL,
            synonyms: self.SYNONYMS_COL,
            inchi_key: self.INCHI_COL
        }
        for param, col in param_to_col_map:
            if param:
                return_set.add(self.translation_map[col].translate(data))

        if len(return_set) != 1:
            return return_set
        return return_set.pop()

    # TODO reverse translate?
