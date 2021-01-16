import os
import csv
import json
from dataclasses import dataclass, field
from mygene import MyGeneInfo


@dataclass()
class DrugBankProteinTarget:
    id: str
    name: str
    gene_name: str
    GenBank_Protein_id: str
    GenBank_gene_id: str
    UniProt_id: str
    UniProt_title: str
    PDB_id: str
    GeneCard_id: str
    GenAtlas_id: str
    HGNC_id: str
    species: str
    drug_ids: list = field(default_factory=list)


class DrugBankProteinTargetsData:
    CSV_ID_KEY = "ID"
    CSV_NAME_KEY = "Name"
    CSV_GENE_NAME_KEY = "Gene Name"
    CSV_GENEBANK_PROT_KEY = "GenBank Protein ID"
    CSV_GENEBANK_GENE_KEY = "GenBank Gene ID"
    CSV_UNIPROT_ID_KEY = "UniProt ID"
    CSV_UNIPROT_TITLE_KEY = "Uniprot Title"
    CSV_PDB_ID_KEY = "PDB ID"
    CSV_GENECARD_ID_KEY = "GeneCard ID"
    CSV_GENEATLAS_ID_KEY = "GenAtlas ID"
    CSV_HGNC_ID_KEY = "HGNC ID"
    CSV_SPECIES_KEY = "Species"
    CSV_DRUG_IDS_KEY = "Drug IDs"

    def __init__(self):
        self._protein_targets = []

    def init_from_file(self, path):
        file_type = os.path.splitext(path)[-1]
        if file_type.lower() == ".csv":
            with open(path, 'r') as handler:
                csv_reader = csv.DictReader(handler)
                self._protein_targets = [
                    DrugBankProteinTarget(id=row[self.CSV_ID_KEY],
                                          name=row[self.CSV_NAME_KEY],
                                          gene_name=row[self.CSV_GENE_NAME_KEY],
                                          GenBank_Protein_id=row[self.CSV_GENEBANK_PROT_KEY],
                                          GenBank_gene_id=row[self.CSV_GENEBANK_GENE_KEY],
                                          UniProt_id=row[self.CSV_UNIPROT_ID_KEY],
                                          UniProt_title=row[self.CSV_UNIPROT_TITLE_KEY],
                                          PDB_id=row[self.CSV_PDB_ID_KEY],
                                          GeneCard_id=row[self.CSV_GENECARD_ID_KEY],
                                          GenAtlas_id=row[self.CSV_GENEATLAS_ID_KEY],
                                          HGNC_id=row[self.CSV_HGNC_ID_KEY],
                                          species=row[self.CSV_SPECIES_KEY],
                                          drug_ids=[d_id.strip() for d_id in row[self.CSV_DRUG_IDS_KEY].split(";")])
                    for row in csv_reader
                ]
        else:
            raise ValueError(f'{self.__class__} can only be initiated from csv files')

    @property
    def protein_targets(self):
        return self._protein_targets

    def get_human_drug_protein_targets(self):
        drug_to_target_map = dict()
        for protein_target in [protein_target for protein_target in self.protein_targets if protein_target.species=="Humans"]:
            protein_symbol = protein_target.gene_name
            targeting_drugs = protein_target.drug_ids
            for drug in targeting_drugs:
                if drug not in drug_to_target_map:
                    drug_to_target_map[drug] = set()
                drug_to_target_map[drug].add(protein_symbol)

        return drug_to_target_map


if __name__ == "__main__":
    with open(r'../../data/gene_symbol_to_entrezgeneid.json', 'r') as json_handler:
        human_gene_symbol_to_id = json.load(json_handler)

    drugbank_data = DrugBankProteinTargetsData()
    drugbank_data.init_from_file(r'C:\studies\thesis\code\NetProp\data\all_drugbank_drug_targets.csv')
    drug_to_target_map = drugbank_data.get_human_drug_protein_targets()
    all_targets = set.union(*[targets for targets in drug_to_target_map.values()])
    ncbi_query = MyGeneInfo().querymany(all_targets, scopes="symbol", fields=["entrezgene", "symbol"], species="human")
    for result in ncbi_query:
        if result.get("notfound", False) or "entrezgene" not in result:
            pass
        else:
            human_gene_symbol_to_id[result["symbol"]] = int(result["entrezgene"])

    with open(r'../../data/gene_symbol_to_entrezgeneid_new.json', 'w') as json_handler:
        human_gene_symbol_to_id = json.dump(human_gene_symbol_to_id, json_handler, indent=4)

