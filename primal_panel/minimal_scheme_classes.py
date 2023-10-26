from Bio import SeqIO
import numpy as np


# Submodule imports
from primal_digest.digestion import digest, generate_valid_primerpairs
from primal_digest.classes import FKmer, RKmer, PrimerPair
from primal_digest.msa import MSA
from primal_digest.mismatches import MatchDB, detect_new_products
from primal_digest.seq_functions import remove_end_insertion, entropy_score_array
from primal_digest.mapping import create_mapping


# Interations checker
from primaldimer_py import do_pools_interact_py

# Other imports
import numpy as np
from abc import ABC, abstractmethod
from uuid import uuid4
from enum import Enum


class Amplicon(ABC):
    @abstractmethod
    def start(self) -> int:
        pass

    @abstractmethod
    def end(self) -> int:
        pass

    @abstractmethod
    def all_seqs(self) -> set[str]:
        pass


class PanelReturn(Enum):
    """
    Enum for the return values of the Panel class.
    """

    ADDED_PRIMERPAIR = 0
    NO_PRIMERPAIRS = 1
    MOVING_TO_NEW_MSA = 2
    NO_MORE_PRIMERPAIRS_IN_MSA = 3


def does_overlap(
    new_pp: tuple[int, int, int], current_pps: list[tuple[int, int, int]]
) -> bool:
    """Checks if a new primerpair overlaps with any of the current primerpair.

    Args:
        new_pp: A tuple representing the new primerpair to check for overlap. The tuple
            contains three integers: the start min(pp.fprimer.starts) , end max(pp.rprimer.ends), and the MSA (multiple sequence
            alignment) index.
        current_pps: A list of tuples representing the current primal panels to check for
            overlap. Each tuple contains three integers: the start index, end index, and the
            MSA index.

    Returns:
        A boolean value indicating whether the new primal panel overlaps with any of the current
        primal panels.
    """
    # Check if from same msa
    for current_pp in current_pps:
        # if not from same msa
        if current_pp[2] != new_pp[2]:
            continue

        # If they overlap
        if range(max(new_pp[0], current_pp[0]), min(new_pp[1], current_pp[1]) + 1):
            return True

    return False


class PanelMSA:
    # Given attributes
    name: str
    path: str
    msa_index: int
    regions: list

    # Calculated on init
    _uuid: str
    _chrom_name: str  # only used in the primer.bed file and html report
    _array: np.ndarray
    _mapping_array: np.ndarray | None

    # Score arrays
    _snp_count_array: np.ndarray | None
    _entropy_array: np.ndarray
    _failed_primerpairs: set[PrimerPair]

    # Calculated attributes
    fkmers: list[FKmer]
    rkmers: list[RKmer]
    primerpairs: list[PrimerPair]

    def __init__(self, name, path, msa_index, mapping, indexes) -> None:
        self.name = name
        self._chrom_name = path.stem
        self.path = path
        self.msa_index = msa_index
        self.indexes = indexes

        # Initialise empty kmer lists
        self.fkmers = []
        self.rkmers = []
        self.primerpairs = []

        # Initialise failed primerpairs
        self._failed_primerpairs = set()

        # Initialise empty arrays
        self.snp_count_array = None

        # Read in the MSA
        records_index = SeqIO.index(str(self.path), "fasta")
        self.array = np.array(
            [record.seq.upper() for record in records_index.values()], dtype="U1"
        )
        self.array = remove_end_insertion(self.array)

        # Create the entropy array
        self._entropy_array = entropy_score_array(self.array)

        # Create the mapping array
        if mapping == "consensus":
            self._mapping_array = None
            self._chrom_name = self.name
        elif mapping == "first":
            self._mapping_array, self.array = create_mapping(self.array, 0)
            self._chrom_name = list(records_index)[0]
        else:
            raise ValueError(f"Unrecognised mapping type: {mapping}")

        # Create the SNP count array if required
        self._snp_count_array = np.zeros(self.array.shape[1], dtype=int)
        if indexes is not None:
            for index in list(indexes):
                self._snp_count_array[index[0] : index[1]] = 1

        # Asign a UUID
        self._uuid = str(uuid4())[:8]

    # override base digest to use indexes
    def digest(self, cfg, indexes=False):
        self.fkmers, self.rkmers = digest(
            msa_array=self.array,
            cfg=cfg,
            indexes=indexes,
        )

    def generate_primerpairs(self, cfg):
        ## Write a doc string
        self.primerpairs = generate_valid_primerpairs(
            self.fkmers,
            self.rkmers,
            cfg,
            self.msa_index,
        )

    def get_pp_entropy(self, pp: PrimerPair) -> float:
        """
        Returns sum of entropy in the primertrimmed amplicon
        """
        return np.sum(self._entropy_array[pp.fprimer.end + 1 : pp.rprimer.start - 1])

    def get_pp_snp_score(self, pp: PrimerPair) -> int:
        """
        Returns number of SNPs in the primertrimmed amplicon
        """
        return np.sum(self._snp_count_array[pp.fprimer.end + 1 : pp.rprimer.start - 1])


class Panel:
    msas: list[MSA]
    cfg: dict
    _pool: list[PrimerPair]
    _matchDB: MatchDB
    _matches: set[tuple]
    _current_msa_index: int

    # for keep adding
    _workingmsasbool: list[bool] | None = None

    def __init__(self, msas: list[MSA], cfg: dict, matchdb: MatchDB) -> None:
        self.msas = msas
        self.cfg = cfg
        self._matchDB = matchdb
        self._matches = set()
        self._pool = []
        self._current_msa_index = 0

    def _next_msa(self) -> int:
        self._current_msa_index = (self._current_msa_index + 1) % len(self.msas)
        return self._current_msa_index

    def _add_primerpair(self, primerpair: PrimerPair, msa_index: int) -> None:
        primerpair.pool = 0
        primerpair.msa_index = msa_index
        primerpair.amplicon_number = len(
            [pp for pp in self._pool if pp.msa_index == msa_index]
        )

        # Adds the primerpair's matches to the pools matches
        self._matches.update(
            primerpair.find_matches(
                self._matchDB,
                fuzzy=self.cfg["mismatch_fuzzy"],
                remove_expected=True,
                kmersize=self.cfg["mismatch_kmersize"],
            )
        )
        # Adds the primerpair to the pool
        self._pool.append(primerpair)

        # Update the current MSA
        self._next_msa()

    def get_covered_regions(self) -> list[tuple[int, int, int]]:
        return [(pp.start, pp.end, pp.msa_index) for pp in self._pool]

    def try_add_primerpair(self) -> bool:
        """
        Try to add a primer pair to the pool.

        :return: True if a primer pair was successfully added to the pool, False otherwise.
        """
        # Get the current MSA
        current_msa = self.msas[self._current_msa_index]

        # Failed PrimerPairs
        failed_primerpairs = set()
        added = False

        all_seqs_in_pool = [
            y for sublist in (x.all_seqs() for x in self._pool) for y in sublist
        ]

        for pos_primerpair in current_msa.primerpairs:
            # Guard if there is overlap
            if does_overlap(
                (pos_primerpair.start, pos_primerpair.end, pos_primerpair.msa_index),
                self.get_covered_regions(),
            ):
                failed_primerpairs.add(pos_primerpair)
                continue

            # Guard if there is an interaction
            if do_pools_interact_py(
                pos_primerpair.all_seqs(),
                all_seqs_in_pool,
                self.cfg["dimerscore"],
            ):
                failed_primerpairs.add(pos_primerpair)
                continue

            # Guard if there is a match
            if detect_new_products(
                pos_primerpair.find_matches(
                    self._matchDB,
                    remove_expected=False,
                    kmersize=self.cfg["mismatch_kmersize"],
                    fuzzy=self.cfg["mismatch_fuzzy"],
                ),
                self._matches,
            ):
                failed_primerpairs.add(pos_primerpair)
                continue

            # If primerpair passes all checks add it to the pool
            self._add_primerpair(pos_primerpair, self._current_msa_index)
            added = True
            break

        # Remove the failed primerpairs
        for failed_primerpair in failed_primerpairs:
            current_msa.primerpairs.remove(failed_primerpair)

        # Return if a primerpair was added
        if added:
            return PanelReturn.ADDED_PRIMERPAIR
        else:
            return PanelReturn.NO_PRIMERPAIRS

    def keep_adding(self):
        # Create a list of which msa still have posible primerpairs
        if self._workingmsasbool is None:
            self._workingmsasbool = [True] * len(self.msas)

        all_seqs_in_pool = [
            y for sublist in (x.all_seqs() for x in self._pool) for y in sublist
        ]

        added = False
        failed_primerpairs: set[PrimerPair] = set()

        if not self._workingmsasbool[self._current_msa_index]:
            self._next_msa()
            return PanelReturn.MOVING_TO_NEW_MSA
        else:
            msa = self.msas[self._current_msa_index]

        # Check each primerpair in this msa
        for pos_primerpair in msa.primerpairs:
            # Guard if there is overlap
            if does_overlap(
                (
                    pos_primerpair.start,
                    pos_primerpair.end,
                    pos_primerpair.msa_index,
                ),
                self.get_covered_regions(),
            ):
                failed_primerpairs.add(pos_primerpair)
                continue

            # Guard if there is an interaction
            if do_pools_interact_py(
                pos_primerpair.all_seqs(),
                all_seqs_in_pool,
                self.cfg["dimerscore"],
            ):
                failed_primerpairs.add(pos_primerpair)
                continue

            # Guard if there is a match
            if detect_new_products(
                pos_primerpair.find_matches(
                    self._matchDB,
                    remove_expected=False,
                    kmersize=self.cfg["mismatch_kmersize"],
                    fuzzy=self.cfg["mismatch_fuzzy"],
                ),
                self._matches,
            ):
                failed_primerpairs.add(pos_primerpair)
                continue

            # If primerpair passes all checks add it to the pool
            self._add_primerpair(pos_primerpair, self._current_msa_index)
            added = True
            break

        # Remove the failed primerpairs
        for failed_primerpair in failed_primerpairs:
            msa.primerpairs.remove(failed_primerpair)

        # Guard for a primerpair was added, return that that primerpair was added
        if added:
            return PanelReturn.ADDED_PRIMERPAIR

        # If no primerpairs were added to this msa, mark it as done
        self._workingmsasbool[self._current_msa_index] = False
        # If there are still working msas, return that we are moving to a new msa
        if any(self._workingmsasbool):
            return PanelReturn.NO_MORE_PRIMERPAIRS_IN_MSA
        else:
            return PanelReturn.NO_PRIMERPAIRS

    def to_bed(self, ref_name: str, amplicon_prefix: str) -> str:
        """
        Returns a bedfile of the current pool
        """
        bed = ""
        for pp in self._pool:
            bed += pp.to_bed(ref_name, amplicon_prefix)
        return bed
