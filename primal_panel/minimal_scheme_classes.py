from Bio import SeqIO
import numpy as np

# Module imports
from primal_panel.entropy import entropy_score_array

# Submodule imports
from primal_digest.digestion import digest, generate_valid_primerpairs
from primal_digest.classes import FKmer, RKmer, PrimerPair
from primal_digest.mismatches import MatchDB, detect_new_products
from primal_digest.seq_functions import remove_end_insertion


# Interations checker
from primaldimer_py import do_pools_interact_py

# Other imports
import itertools
import numpy as np
from abc import ABC, abstractmethod, abstractproperty
import random


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


class NestedFKmer:
    inner: FKmer
    outer: FKmer

    start: int
    end: int

    def __init__(self, inner: FKmer, outer: FKmer) -> None:
        self.inner = inner
        self.outer = outer

        self.start = min(outer.starts())
        self.end = inner.end

    def all_seqs(self) -> set[str]:
        return self.inner.seqs | self.outer.seqs


class NestedRKmer:
    inner: RKmer
    outer: RKmer

    start: int
    end: int

    def __init__(self, inner: RKmer, outer: RKmer) -> None:
        self.inner = inner
        self.outer = outer

        self.start = outer.start
        self.end = max(inner.ends())

    def all_seqs(self) -> set[str]:
        return self.inner.seqs | self.outer.seqs


class NestedPP(Amplicon):
    outerfkmer: FKmer
    outerrkmer: RKmer
    primerpair: PrimerPair

    def __init__(self, fkmer, rkmer, primerpair) -> None:
        self.outerfkmer = fkmer
        self.outerrkmer = rkmer
        self.primerpair = primerpair

    def all_seqs(self) -> set[str]:
        return (
            self.outerfkmer.seqs
            | self.outerrkmer.seqs
            | set(self.primerpair.all_seqs())
        )

    def start(self) -> int:
        return self.outerfkmer.start

    def end(self) -> int:
        return self.outerrkmer.end

    def __str__(self, ref_name, amplicon_prefix):
        return (
            self.outerfkmer.__str__(
                referance=f"{ref_name}",
                amplicon_prefix=f"{amplicon_prefix}_OUTER_{0}",
                pool=0,
            )
            + self.primerpair.__str__(ref_name, f"{amplicon_prefix}_INNER")
            + self.outerrkmer.__str__(
                referance=f"{ref_name}",
                amplicon_prefix=f"{amplicon_prefix}_OUTER_{0}",
                pool=0,
            )
        )


def does_overlap(
    new_pp: tuple[int, int, int], current_pps: list[tuple[int, int, int]]
) -> bool:
    """
    Returns true if the Amplicons overlap
    Input should be a list of tuples of (start, end, msa_index)
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

        # Create the entropy array
        self._entropy_array = entropy_score_array(self.array)

    def digest(self, cfg, indexes=False):
        self.fkmers, self.rkmers = digest(
            msa_array=self.array,
            cfg=cfg,
            indexes=indexes,
        )

    def generate_primerpairs(self, cfg):
        self.primerpairs = generate_valid_primerpairs(
            self.fkmers,
            self.rkmers,
            cfg,
            self.msa_index,
        )

    def get_pp_entropy(self, pp: PrimerPair) -> float:
        return np.sum(self._entropy_array[pp.fprimer.end + 1 : pp.rprimer.start - 1])


class panel:
    msas: list[MSA]
    cfg: dict
    _pool: list[PrimerPair]
    _matchDB: MatchDB
    _matches: set[tuple]
    _current_msa_index: int

    def __init__(self, msas: list[MSA], cfg: dict, matchdb: MatchDB) -> None:
        self.msas = msas
        self.cfg = cfg
        self._matchDB = matchdb
        self._matches = set()
        self._pool = []
        self._current_msa_index = 0

    def _next_msa(self) -> None:
        self._current_msa_index = (self._current_msa_index + 1) % len(self.msas)

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

        all_seqs_in_pool = [
            y for sublist in (x.all_seqs() for x in self._pool) for y in sublist
        ]

        for pos_primerpair in current_msa.primerpairs:
            # Guard if there is overlap
            if does_overlap(
                (pos_primerpair.start, pos_primerpair.end, pos_primerpair.msa_index),
                self.get_covered_regions(),
            ):
                continue

            # Guard if there is an interaction
            if do_pools_interact_py(
                pos_primerpair.all_seqs(),
                all_seqs_in_pool,
                self.cfg["dimerscore"],
            ):
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
                continue

            # If primerpair passes all checks add it to the pool
            self._add_primerpair(pos_primerpair, self._current_msa_index)
            return True

        # If no primerpairs were added
        return False

    def to_bed(self, ref_name: str, amplicon_prefix: str) -> str:
        """
        Returns a bedfile of the current pool
        """
        bed = ""
        for pp in self._pool:
            bed += pp.to_bed(ref_name, amplicon_prefix)
        return bed
