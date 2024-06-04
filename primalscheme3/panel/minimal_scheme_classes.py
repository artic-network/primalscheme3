import pathlib
from enum import Enum
from itertools import islice

import numpy as np

# Interations checker
from primaldimer_py import do_pools_interact_py  # type: ignore

from primalscheme3.core.bedfiles import BedPrimerPair

# Core Module imports
from primalscheme3.core.classes import PrimerPair
from primalscheme3.core.mismatches import MatchDB, detect_new_products
from primalscheme3.core.msa import MSA
from primalscheme3.core.multiplex import Multiplex
from primalscheme3.core.seq_functions import (
    entropy_score_array,
)
from primalscheme3.core.progress_tracker import ProgressManager


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


class Region:
    chromname: str
    start: int
    stop: int
    name: str
    score: int

    def __init__(self, chromanem: str, start: int, stop: int, name: str, score: int):
        self.chromname = chromanem
        self.start = start
        self.stop = stop
        self.name = name
        self.score = score

    def positions(self):
        return range(self.start, self.stop)

    def __hash__(self) -> int:
        return hash(
            f"{self.chromname}:{self.start}:{self.stop}:{self.name}:{self.score}"
        )

    def __eq__(self, __value: object) -> bool:
        if not isinstance(__value, Region):
            return False
        return hash(self) == hash(__value)

    def to_bed(self) -> str:
        return f"{self.chromname}\t{self.start}\t{self.stop}\t{self.name}\t{self.score}"


class PanelMSA(MSA):
    # Score arrays
    _score_array: np.ndarray
    _entropy_array: np.ndarray

    primerpairpointer: int

    def __init__(
        self,
        name: str,
        path: pathlib.Path,
        msa_index: int,
        mapping: str,
        progress_manager: ProgressManager,
        logger=None,
    ) -> None:
        # Call the MSA init
        super().__init__(
            name=name,
            path=path,
            msa_index=msa_index,
            mapping=mapping,
            logger=logger,
            progress_manager=progress_manager,
        )

        # Create the entropy array
        self._entropy_array = np.array(entropy_score_array(self.array))

        # Create the ref_to_msa_index_array
        self.ref_to_msa_index_array: list[None | int] = [None] * self.array.shape[1]
        if self._mapping_array is not None:
            for msa_index, ref_index in enumerate(self._mapping_array):
                if ref_index is not None:
                    self.ref_to_msa_index_array[ref_index] = msa_index

        # Create the primerpairpointer
        self.primerpairpointer = 0
        self.regions = None

    def add_regions(self, regions: list[Region]) -> None:
        self.regions = regions

        # Create the score array
        self._score_array = np.zeros(self.array.shape[1], dtype=int)
        for region in regions:
            # Map the ref positions to msa positions
            msa_region_start = self.ref_to_msa_index_array[region.start]
            msa_region_stop = self.ref_to_msa_index_array[region.stop]
            # Add the score to the score array
            if msa_region_start is not None and msa_region_stop is not None:
                for i in range(msa_region_start, msa_region_stop):
                    self._score_array[i] += region.score

    def remove_kmers_that_clash_with_regions(self):
        """
        Removes f/rkmers who clash with the regions
        """
        if self.regions is None:
            return

        # Remove primer that overlap with regions
        regions = [(x.start, x.stop, self.msa_index) for x in self.regions]
        self.fkmers = [
            fkmer
            for fkmer in self.fkmers
            if not does_overlap(
                (min(fkmer.starts()), fkmer.end, self.msa_index),
                regions,
            )
        ]
        self.rkmers = [
            rkmer
            for rkmer in self.rkmers
            if not does_overlap(
                (rkmer.start, max(rkmer.ends()), self.msa_index),
                regions,
            )
        ]

    def iter_unchecked_primerpairs(self):
        """
        Returns all primerpairs that have not been checked yet
        """
        return islice(self.primerpairs, self.primerpairpointer, None)

    def get_pp_entropy(self, pp: PrimerPair) -> float:
        """
        Returns sum of entropy in the primertrimmed amplicon
        """
        return np.sum(self._entropy_array[pp.fprimer.end : pp.rprimer.start - 1])

    def get_pp_score(self, pp: PrimerPair) -> int:
        """
        Returns number of SNPs in the primertrimmed amplicon
        """
        msa_pp_fprimer_end = self.ref_to_msa_index_array[pp.fprimer.end]
        msa_pp_rprimer_start = self.ref_to_msa_index_array[pp.rprimer.start]
        if msa_pp_fprimer_end is None or msa_pp_rprimer_start is None:
            return None
        return np.sum(self._score_array[msa_pp_fprimer_end : msa_pp_rprimer_start - 1])


class Panel(Multiplex):
    # Base class
    _pools: list[list[PrimerPair | BedPrimerPair]]
    _current_pool: int
    _last_pp_added: list[PrimerPair]  # Stack to keep track of the last primer added
    _matchDB: MatchDB
    cfg: dict

    # New attributes
    msas: list[PanelMSA]
    _current_msa_index: int

    # for keep adding
    _workingmsasbool: list[bool] | None = None

    def __init__(self, msas: list[PanelMSA], cfg: dict, matchdb: MatchDB) -> None:
        super().__init__(cfg, matchdb)
        self.msas = msas
        self._current_msa_index = 0
        self._failed_primerpairs = [set() for _ in range(self.n_pools)]

    def _next_msa(self) -> int:
        """
        Updates the current msa index to the next msa. Returns the new msa index.
        :return: The new msa index.
        """
        self._current_msa_index = (self._current_msa_index + 1) % len(self.msas)
        return self._current_msa_index

    def _add_primerpair(
        self, primerpair: PrimerPair, pool: int, msa_index: int
    ) -> None:
        # Add the primerpair to the pool
        super().add_primer_pair_to_pool(primerpair, pool, msa_index)
        # Update the current MSA
        self._next_msa()

    def get_covered_regions(self, pool: int) -> list[tuple[int, int, int]]:
        """
        Returns all regions covered in the pool
        """
        return [(pp.start, pp.end, pp.msa_index) for pp in self._pools[pool]]

    def try_add_primerpair(self) -> PanelReturn:
        """
        Try to add a primer pair to the current pool.
        :return: PanelReturn.
        """
        # Get the current MSA
        current_msa = self.msas[self._current_msa_index]
        current_pool = self._current_pool

        # Get the current pool
        added = False
        # All seqs in each pool
        seqs_in_pool = [
            seq
            for seq in (pp.all_seqs() for pp in self._pools[current_pool])
            for seq in seq
        ]
        # For each primerpair in the current msa
        for new_pointer, pos_primerpair in enumerate(
            current_msa.iter_unchecked_primerpairs(), current_msa.primerpairpointer
        ):
            # Guard if the primerpair is in the failed primerpairs
            if pos_primerpair in self._failed_primerpairs[current_pool]:
                continue

            # Guard if there is overlap
            if super().does_overlap(pos_primerpair, current_pool):
                self._failed_primerpairs[current_pool].add(pos_primerpair)
                continue

            # Guard if there is an interaction
            if do_pools_interact_py(
                pos_primerpair.all_seqs(),
                seqs_in_pool,
                self.cfg["dimerscore"],
            ):
                self._failed_primerpairs[current_pool].add(pos_primerpair)
                continue

            # Guard if there is a match
            if detect_new_products(
                pos_primerpair.find_matches(
                    self._matchDB,
                    remove_expected=False,
                    kmersize=self.cfg["mismatch_kmersize"],
                    fuzzy=self.cfg["mismatch_fuzzy"],
                ),
                self._matches[current_pool],
            ):
                self._failed_primerpairs[current_pool].add(pos_primerpair)
                continue

            # If primerpair passes all checks add it to the pool
            self._add_primerpair(pos_primerpair, current_pool, self._current_msa_index)
            added = True
            # Update the pointer
            current_msa.primerpairpointer = new_pointer
            break

        # Return if a primerpair was added
        if added:
            return PanelReturn.ADDED_PRIMERPAIR
        else:
            return PanelReturn.NO_PRIMERPAIRS

    def keep_adding(self) -> PanelReturn:
        # Create a list of which msa still have posible primerpairs
        if self._workingmsasbool is None:
            self._workingmsasbool = [True] * len(self.msas)

        # All seqs in each pool
        all_seqs_in_pools: dict[int, list[str]] = {
            pool: [
                y
                for sublist in (x.all_seqs() for x in self._pools[pool])
                for y in sublist
            ]
            for pool in range(self.n_pools)
        }
        # Move to next msa if current msa is done
        if not self._workingmsasbool[self._current_msa_index]:
            self._next_msa()
            return PanelReturn.MOVING_TO_NEW_MSA
        else:
            msa = self.msas[self._current_msa_index]

        # For each primerpair
        for new_pointer, pos_primerpair in enumerate(
            msa.iter_unchecked_primerpairs(), msa.primerpairpointer
        ):
            # For each pool
            for pool in range(self.n_pools):
                # Guard if there is overlap
                if super().does_overlap(pos_primerpair, pool):
                    continue

                # Guard if there is an interaction
                if do_pools_interact_py(
                    pos_primerpair.all_seqs(),
                    all_seqs_in_pools[pool],
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
                    self._matches[pool],
                ):
                    continue

                # If primerpair passes all checks add it to the pool
                self._add_primerpair(pos_primerpair, pool, self._current_msa_index)
                # Update the pointer
                msa.primerpairpointer = new_pointer
                return PanelReturn.ADDED_PRIMERPAIR

        # If no primerpairs were added to this msa, mark it as done
        self._workingmsasbool[self._current_msa_index] = False
        # If there are still working msas, return that we are moving to a new msa
        if any(self._workingmsasbool):
            return PanelReturn.NO_MORE_PRIMERPAIRS_IN_MSA
        else:
            return PanelReturn.NO_PRIMERPAIRS
