import numpy as np
from Bio import SeqIO


# Module imports
from primal_digest.classes import FKmer, RKmer, PrimerPair
from primal_digest.seq_functions import remove_end_insertion
from primal_digest.digestion import digest, generate_valid_primerpairs


class MSA:
    name: str
    path: str
    msa_index: int
    array: np.ndarray
    fkmers: list[FKmer]
    rkmers: list[RKmer]
    primerpairs: list[PrimerPair]

    def __init__(self, name, path, msa_index) -> None:
        self.name = name
        self.path = path
        self.msa_index = msa_index

        # Read in the MSA
        records = SeqIO.parse(self.path, "fasta")
        self.array = np.array([record.seq.upper() for record in records], dtype=str)
        self.array = remove_end_insertion(self.array)

    def digest(self, cfg):
        self.fkmers, self.rkmers = digest(self.array, cfg)

    def generate_primerpairs(self, cfg):
        self.primerpairs = generate_valid_primerpairs(
            self.fkmers, self.rkmers, cfg, self
        )
