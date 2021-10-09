from dataclasses import dataclass, field
from propagater import PropagationContainer
from pydantic import BaseModel
from typing import List



HUMAN_SPECIES_NAME = "human"
COVID_19_SPECIES_NAME = "covid"

@dataclass()
class Protein(PropagationContainer):
    species: str = HUMAN_SPECIES_NAME


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


@dataclass(frozen=True)
class GeneSet:
    name: str
    gene_set: set = field(default_factory=set)

    def __iter__(self):
        return iter(self.gene_set)


class LiquidModel(BaseModel):
    sources: List
    targets: List


