import os
import pandas as pd

class DrugbankData:
    def __init__(self):
        self.drugbank_dataframe = None

    def init_from_file(self, path):
        file_type = os.path.splitext(path)[-1]
        if file_type.lower() == ".csv":
            self.drugbank_dataframe = pd.read_csv(path, index_col=None)
        else:
            raise ValueError(f'{self.__class__} can only be initiated from csv files')


class ProteinTargetData(DrugbankData):
    _PROTEIN_ID_COL = "GenBank Protein ID"
    _DRUG_IDS_COL = "Drug IDs"
    _DRUG_NAME_SEPERATOR = ";"

    def __init__(self):
        super().__init__()
        self.drug_targets = dict()

    def init_from_file(self, path):
        super().init_from_file(path)
        self._drug_targets_from_dataframe()

    def _drug_targets_from_dataframe(self):
        for row in self.drugbank_dataframe.iterrows():
            protein_id = row[1][self._PROTEIN_ID_COL]
            for drug_name in row[1][self._DRUG_IDS_COL].split(self._DRUG_NAME_SEPERATOR):
                if drug_name not in self.drug_targets:
                    self.drug_targets[drug_name] = {str(protein_id)}
                else:
                    self.drug_targets[drug_name].add(str(protein_id))

d = ProteinTargetData()
d.init_from_file(r'D:\Etay\studies\thesis\drugbank_data\all.csv')

