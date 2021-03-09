from dataclasses import dataclass, field
from functools import total_ordering


@dataclass(frozen=False)
@total_ordering
class Protein:
    id: int
    source_of: set = field(default_factory=set)
    target_of: set = field(default_factory=set)

    def __hash__(self):
        return self.id

    def __eq__(self, other):
        return self.id == other.id

    def __lt__(self, other):
        return self.liquid < other.liquid


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