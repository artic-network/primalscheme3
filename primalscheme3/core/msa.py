import numpy as np
from Bio import SeqIO
from uuid import uuid4

from primalscheme3.core.digestion import digest, generate_valid_primerpairs
from primalscheme3.core.seq_functions import remove_end_insertion
from primalscheme3.core.classes import FKmer, RKmer, PrimerPair
from primalscheme3.core.mapping import create_mapping


class MSA:
    # Provided
    name: str
    path: str
    msa_index: int

    # Calculated on init
    array: np.ndarray
    _uuid: str
    _chrom_name: str  # only used in the primer.bed file and html report
    _mapping_array: np.ndarray | None

    # Calculated on evaluation
    fkmers: list[FKmer]
    rkmers: list[RKmer]
    primerpairs: list[PrimerPair]

    def __init__(self, name, path, msa_index, mapping) -> None:
        self.name = name
        self.path = str(path)
        self.msa_index = msa_index

        # Read in the MSA
        records_index = SeqIO.index(self.path, "fasta")
        self.array = np.array(
            [record.seq.upper() for record in records_index.values()], dtype="U1"
        )
        self.array = remove_end_insertion(self.array)

        # Create the mapping array
        if mapping == "consensus":
            self._mapping_array = None
            self._chrom_name = self.name
        elif mapping == "first":
            self._mapping_array, self.array = create_mapping(self.array, 0)
            self._chrom_name = list(records_index)[0]

        # Asign a UUID
        self._uuid = str(uuid4())[:8]

    def digest(
        self, cfg: dict, indexes: tuple[list[int], list[int]] | None = None
    ) -> None:
        """
        Digest the given MSA array and return the FKmers and RKmers.

        :param cfg: A dictionary containing configuration parameters.
        :param indexes: A tuple of MSA indexes for (FKmers, RKmers), or False to use all indexes.
        :return: None (Class is updated inplace)
        """
        # Create all the kmers
        self.fkmers, self.rkmers = digest(
            msa_array=self.array,
            cfg=cfg,
            indexes=indexes,
        )
        # remap the fkmer and rkmers if needed
        if self._mapping_array is not None:
            self.fkmers = [fkmer.remap(self._mapping_array) for fkmer in self.fkmers]  # type: ignore
            self.fkmers = [x for x in self.fkmers if x is not None]
            self.rkmers = [rkmer.remap(self._mapping_array) for rkmer in self.rkmers]  # type: ignore
            self.rkmers = [x for x in self.rkmers if x is not None]

    def generate_primerpairs(self, cfg) -> None:
        self.primerpairs = generate_valid_primerpairs(
            self.fkmers,
            self.rkmers,
            cfg,
            self.msa_index,
        )
