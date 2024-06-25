import pathlib
import sys
from uuid import uuid4

import numpy as np
from Bio import SeqIO

from primalscheme3.core.classes import FKmer, PrimerPair, RKmer
from primalscheme3.core.digestion import digest, generate_valid_primerpairs
from primalscheme3.core.mapping import create_mapping
from primalscheme3.core.seq_functions import remove_end_insertion


def msa_qc(array: np.ndarray):
    """
    Checks all sequences are the same length. Checks there are no empty cols
    """
    empty_set = {"", "-"}

    # Check rows are same length
    for row_index in range(0, array.shape[0]):
        if len(array[row_index, :]) != len(array[0]):
            raise ValueError("MSA contains sequences of different lengths")

    # Check for empty columns
    for col_index in range(0, array.shape[1]):
        slice: set[str] = set(array[:, col_index])
        if slice.issubset(empty_set):
            raise ValueError(f"MSA contains empty column at: {col_index}")


class MSA:
    # Provided
    name: str
    path: str
    msa_index: int

    # Calculated on init
    array: np.ndarray
    _uuid: str
    _chrom_name: str  # only used in the primer.bed file and html report
    _mapping_array: np.ndarray
    _ref_to_msa: dict[int, int]

    # Calculated on evaluation
    fkmers: list[FKmer]
    rkmers: list[RKmer]
    primerpairs: list[PrimerPair]

    def __init__(
        self,
        name: str,
        path: pathlib.Path,
        msa_index: int,
        mapping: str,
        progress_manager,
        logger=None,
    ) -> None:
        self.name = name
        self.path = str(path)
        self.msa_index = msa_index
        self.logger = logger
        self.progress_manager = progress_manager

        # Read in the MSA
        records_index = SeqIO.index(self.path, "fasta")
        self.array = np.array(
            [record.seq.upper() for record in records_index.values()],
            dtype="U1",
            ndmin=2,  # Enforce 2D array even if one genome
        )
        # Do some basic QC
        try:
            msa_qc(self.array)
        except ValueError as e:
            if self.logger:
                self.logger.error(f"MSA: {self.name} failed QC: {e}")
            else:
                print(f"MSA: {self.name} failed QC: {e}")
            sys.exit(1)

        self.array = remove_end_insertion(self.array)

        # Create the mapping array
        # Goes from msa idx -> ref idx
        if mapping == "consensus":
            self._chrom_name = self.name + "_consensus"
            self._mapping_array = np.array([*range(len(self.array[0]))])
        elif mapping == "first":
            self._mapping_array, self.array = create_mapping(self.array, 0)
            self._chrom_name = list(records_index)[0]
        else:
            raise ValueError(f"Mapping method: {mapping} not recognised")

        # Create the Unmapping array
        # Goes from ref idx -> msa idx
        self._ref_to_msa = {
            x: i for i, x in enumerate(list(self._mapping_array)) if x is not None
        }

        # Asign a UUID
        self._uuid = str(uuid4())[:8]

    def digest(
        self,
        cfg: dict,
        indexes: tuple[list[int], list[int]] | None = None,
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
            logger=self.logger,
            progress_manager=self.progress_manager,
            chrom=self.name,
        )
        # remap the fkmer and rkmers if needed
        if self._mapping_array is not None:
            mapping_set = set(self._mapping_array)

            remaped_fkmers = [fkmer.remap(self._mapping_array) for fkmer in self.fkmers]  # type: ignore
            self.fkmers = [
                x
                for x in remaped_fkmers
                if x is not None and x.end in mapping_set and min(x.starts()) >= 0
            ]
            remaped_rkmers = [rkmer.remap(self._mapping_array) for rkmer in self.rkmers]  # type: ignore
            self.rkmers = [
                x
                for x in remaped_rkmers
                if x is not None
                and x.start in mapping_set
                and max(x.ends()) < self.array.shape[1]
            ]

    def generate_primerpairs(
        self, amplicon_size_min: int, amplicon_size_max: int, dimerscore: float
    ) -> None:
        self.primerpairs = generate_valid_primerpairs(
            fkmers=self.fkmers,
            rkmers=self.rkmers,
            amplicon_size_min=amplicon_size_min,
            amplicon_size_max=amplicon_size_max,
            dimerscore=dimerscore,
            msa_index=self.msa_index,
            progress_manager=self.progress_manager,
            chrom=self.name,
        )
        # Update primerpairs to include the chrom_name and amplicon_prefix
        for primerpair in self.primerpairs:
            primerpair.chrom_name = self._chrom_name
            primerpair.amplicon_prefix = self._uuid
