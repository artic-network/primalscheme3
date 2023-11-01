# Module imports
from primalscheme3.core.seq_functions import reverse_complement, remove_end_insertion
from primalscheme3.core.mismatches import MatchDB
from primalscheme3.core.digestion import (
    digest,
    generate_valid_primerpairs,
)
from primalscheme3.core.mapping import create_mapping

from primaldimer_py import do_pools_interact_py


class FKmer:
    end: int
    seqs: set[str]
    _starts: set[int]

    # Add slots for some performance gains
    __slots__ = ["end", "seqs", "_starts"]

    def __init__(self, end, seqs) -> None:
        self.end = end
        self.seqs = seqs
        self._starts = {self.end - len(x) for x in self.seqs}

    def len(self) -> set[int]:
        return {len(x) for x in self.seqs}

    def starts(self) -> set[int]:
        return self._starts

    def __str__(self, reference, amplicon_prefix, pool) -> str:
        string_list = []
        counter = 0
        seqs = list(self.seqs)
        seqs.sort()
        for seq in seqs:
            string_list.append(
                f"{reference}\t{self.end-len(seq)}\t{self.end}\t{amplicon_prefix}_LEFT_{counter}\t{pool}\t+\t{seq}\n"
            )
            counter += 1
        return "".join(string_list)

    def find_matches(
        self,
        matchDB: MatchDB,
        remove_expected: bool,
        fuzzy: bool,
        kmersize: int,
        msa_index,
    ) -> set[tuple]:
        """Returns all matches of this FKmer"""
        return matchDB.find_fkmer(
            self,
            fuzzy=fuzzy,
            remove_expected=remove_expected,
            kmersize=kmersize,
            msaindex=msa_index,
        )

    def __hash__(self) -> int:
        seqs = list(self.seqs)
        seqs.sort()
        return hash(f"{self.end}{self.seqs}")

    def __eq__(self, other):
        if isinstance(other, FKmer):
            return self.__hash__() == other.__hash__()
        else:
            return False

    def remap(self, mapping_array):
        """
        Remaps the fkmer to a new indexing system
        Returns None if the fkmer is not valid
        """
        if mapping_array[self.end] is not None:
            self.end = mapping_array[self.end]
            self._starts = {self.end - len(x) for x in self.seqs}
            return self
        else:
            return None


class RKmer:
    start: int
    seqs: set[str]
    _ends: set[int]

    # Add slots for some performance gains
    __slots__ = ["start", "seqs", "_ends"]

    def __init__(self, start, seqs) -> None:
        self.start = start
        self.seqs = seqs
        self._ends = {len(x) + self.start for x in self.seqs}

    def len(self) -> set[int]:
        return {len(x) for x in self.seqs}

    def ends(self) -> set[int]:
        return self._ends

    def __str__(self, reference, amplicon_prefix, pool) -> str:
        string_list = []
        counter = 0
        seqs = list(self.seqs)
        seqs.sort()
        for seq in seqs:
            string_list.append(
                f"{reference}\t{self.start}\t{self.start+len(seq)}\t{amplicon_prefix}_RIGHT_{counter}\t{pool}\t-\t{seq}\n"
            )
            counter += 1
        return "".join(string_list)

    def reverse_complement(self) -> set[str]:
        return {reverse_complement(x) for x in self.seqs}

    def find_matches(
        self,
        matchDB: MatchDB,
        remove_expected: bool,
        fuzzy: bool,
        kmersize: int,
        msa_index: int,
    ) -> set[tuple]:
        """Returns all matches of this FKmer"""
        return matchDB.find_rkmer(
            self,
            fuzzy=fuzzy,
            remove_expected=remove_expected,
            kmersize=kmersize,
            msaindex=msa_index,
        )

    def __hash__(self) -> int:
        seqs = list(self.seqs)
        seqs.sort()
        return hash(f"{self.start}{self.seqs}")

    def __eq__(self, other):
        if isinstance(other, RKmer):
            return self.__hash__() == other.__hash__()
        else:
            return False

    def remap(self, mapping_array):
        """
        Remaps the rkmer to a new indexing system
        Returns None if the rkmer is not valid
        """
        if mapping_array[self.start] is not None:
            self.start = mapping_array[self.start]
            self._ends = {len(x) + self.start for x in self.seqs}
            return self
        else:
            return None


class PrimerPair:
    fprimer: FKmer
    rprimer: RKmer
    amplicon_number: int
    pool: int
    msa_index: int

    __slots__ = ["fprimer", "rprimer", "amplicon_number", "pool", "msa_index"]

    def __init__(
        self,
        fprimer,
        rprimer,
        msa_index,
        amplicon_number=-1,
        pool=-1,
    ):
        self.fprimer = fprimer
        self.rprimer = rprimer
        self.amplicon_number = amplicon_number
        self.pool = pool
        self.msa_index = msa_index

    def set_amplicon_number(self, amplicon_number) -> None:
        self.amplicon_number = amplicon_number

    def set_pool_number(self, pool_number) -> None:
        self.amplicon_number = pool_number

    def find_matches(self, matchDB, fuzzy, remove_expected, kmersize) -> set[tuple]:
        """
        Find matches for the FKmer and RKmer
        """
        matches = set()
        # Find the FKmer matches
        matches.update(
            self.fprimer.find_matches(
                matchDB, fuzzy, remove_expected, kmersize, msa_index=self.msa_index
            )
        )
        # Find the RKmer matches
        matches.update(
            self.rprimer.find_matches(
                matchDB, fuzzy, remove_expected, kmersize, self.msa_index
            )
        )
        return matches

    @property
    def start(self) -> int:
        return min(self.fprimer.starts())

    @property
    def end(self) -> int:
        return max(self.rprimer.ends())

    def inter_free(self, cfg) -> bool:
        """
        True means interaction
        """
        return do_pools_interact_py(
            self.fprimer.seqs, self.rprimer.seqs, cfg["dimerscore"]
        )

    def all_seqs(self) -> set[str]:
        return [x for x in self.fprimer.seqs] + [x for x in self.rprimer.seqs]

    def __hash__(self) -> int:
        return hash(f"{self.start}{self.end}{self.all_seqs()}")

    def __str__(self, ref_name, amplicon_prefix):
        return self.fprimer.__str__(
            reference=f"{ref_name}",
            amplicon_prefix=f"{amplicon_prefix}_{self.amplicon_number}",
            pool=self.pool + 1,
        ) + self.rprimer.__str__(
            reference=f"{ref_name}",
            amplicon_prefix=f"{amplicon_prefix}_{self.amplicon_number}",
            pool=self.pool + 1,
        )


import numpy as np
from Bio import SeqIO
from uuid import uuid4


class MSA:
    # Provided
    name: str
    path: str
    msa_index: int
    array: np.ndarray

    # Calculated on init
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

    def digest(self, cfg, indexes: bool = False):
        self.fkmers, self.rkmers = digest(
            msa_array=self.array,
            cfg=cfg,
            indexes=indexes,
        )
        # remap the fkmer and rkmers if needed
        if self._mapping_array is not None:
            self.fkmers = [fkmer.remap(self._mapping_array) for fkmer in self.fkmers]
            self.fkmers = [x for x in self.fkmers if x is not None]
            self.rkmers = [rkmer.remap(self._mapping_array) for rkmer in self.rkmers]
            self.rkmers = [x for x in self.rkmers if x is not None]

    def generate_primerpairs(self, cfg):
        self.primerpairs = generate_valid_primerpairs(
            self.fkmers,
            self.rkmers,
            cfg,
            self.msa_index,
        )
