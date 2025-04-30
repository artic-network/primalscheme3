from enum import Enum

import numpy as np
from primalschemers._core import do_pool_interact  # type: ignore

from primalscheme3.core.bedfiles import (
    BedPrimerPair,
    create_amplicon_str,
    create_bedfile_str,
)
from primalscheme3.core.classes import PrimerPair
from primalscheme3.core.config import Config
from primalscheme3.core.mismatches import MatchDB, detect_new_products
from primalscheme3.core.msa import MSA


class PrimerPairCheck(Enum):
    """
    Enum for checking if a primerpair can be added to a pool
    """

    OVERLAP = 1
    INTERACTING = 2
    MISPRIMING = 3
    OK = 4
    CIRCULAR = 5


class Multiplex:
    """
    This is the baseclass for all multiplexes (Scheme / Panel)
    - It allows multiple pools

    Params:
    - cfg: dict. The configuration dict
    - matchDB: MatchDB object
    - msa_dict: dict[int, MSA]. A dict of MSA objects

    Internal:
    - self._pools: list[list[PrimerPair | BedPrimerPair]]. List of pools with PrimerPairs
    - self._current_pool: int. The current pool number
    - self._last_pp_added: list[PrimerPair]. Stack to keep track of the last primer added
    - self._matchDB: MatchDB object
    - self._matches: list[set[tuple]]. A list of sets with matches for each pool
    - self._msa_dict: dict[int, MSA]. A dict of MSA objects
    - self._coverage: dict[int, np.ndarray]. PrimerTrimmed Coverage for each MSA index
    - self._lookup: dict[int, np.ndarray[npools, ncols]] | None. A lookup for the primerpairs in the multiplex.
    """

    _pools: list[list[PrimerPair | BedPrimerPair]]
    _current_pool: int
    _last_pp_added: list[PrimerPair]  # Stack to keep track of the last primer added
    _matchDB: MatchDB
    _matches: list[set[tuple]]
    _msa_dict: dict[int, MSA]
    _coverage: dict[int, np.ndarray]  # PrimerTrimmed Coverage for each MSA index
    _lookup: dict[int, np.ndarray]
    config: Config

    def __init__(
        self, config: Config, matchDB: MatchDB, msa_dict: dict[int, MSA]
    ) -> None:
        self.n_pools = config.n_pools
        self._pools = [[] for _ in range(self.n_pools)]
        self._matches: list[set[tuple]] = [set() for _ in range(self.n_pools)]
        self._current_pool = 0
        self._pp_number = 1
        self.config = config
        self._matchDB = matchDB
        self._last_pp_added = []
        self._msa_dict = msa_dict
        # Set up coverage dict
        self.setup_coverage()

        # Set up the lookup
        self.setup_primerpair_lookup()

    def setup_primerpair_lookup(self) -> None:
        """
        Returns a lookup of primerpairs in the multiplex
        :return: None
        """
        self._lookup = {}
        for msa_index, msa in self._msa_dict.items():
            if msa._mapping_array is None:
                n = len(msa.array[1])
            else:
                n = len(msa._mapping_array)

            n = [None] * n
            # Create a lookup for the primerpairs
            self._lookup[msa_index] = np.array(
                [n for _ in range(self.n_pools)], ndmin=2
            )

    def update_lookup(self, primerpair: PrimerPair | BedPrimerPair, add: bool) -> None:
        """
        Updates the lookup for the primerpair
        :param primerpair: PrimerPair object
        :param add: bool. If True, add to the lookup. If False, remove from the lookup
        :return: None
        """
        circular = primerpair.fprimer.end > primerpair.rprimer.start

        # If the msa_index is not in the lookup. Then return
        if primerpair.msa_index not in self._lookup:
            return

        value = primerpair if add else None

        if circular:
            self._lookup[primerpair.msa_index][
                primerpair.pool, primerpair.fprimer.region()[0] :
            ] = value
            self._lookup[primerpair.msa_index][
                primerpair.pool, : primerpair.rprimer.region()[1]  # slice non inclusive
            ] = value
        else:
            self._lookup[primerpair.msa_index][
                primerpair.pool,
                primerpair.fprimer.region()[0] : primerpair.rprimer.region()[
                    1
                ],  # slice non inclusive
            ] = value

    def get_seqs_in_pool(self, pool: int) -> list[str]:
        """
        Returns a list of all the sequences in the pool
        :param pool: int
        :return: list[str]
        """
        return [
            seq
            for ppseq in (pp.all_seqs() for pp in self._pools[pool])
            for seq in ppseq
        ]

    def get_seqs_bytes_in_pool(self, pool: int) -> list[bytes]:
        """
        Returns a list of all the sequences in the pool
        :param pool: int
        :return: list[str]
        """
        return [
            seq
            for ppseq in (pp.all_seq_bytes() for pp in self._pools[pool])
            for seq in ppseq
        ]

    def check_primerpair_can_be_added(
        self,
        primerpair: PrimerPair,
        pool: int,
        otherseqs_bytes: list[bytes] | None = None,
    ) -> PrimerPairCheck:
        """
        Checks if the primerpair can be added to the multiplex
        - Checks if the primerpair overlaps with any other primerpairs in the pool
        - Checks if the primerpair has any matches in the MatchDB
        :param primerpair: PrimerPair object
        :param otherseqs: list[str] | None. A list of other sequences to check for interactions against
        :return: bool. True if the primerpair can be added
        """
        # Check if the primerpair overlaps with any other primerpairs in the pool
        if self.does_overlap(primerpair, pool):
            return PrimerPairCheck.OVERLAP

        # Check for interactions with other sequences
        if otherseqs_bytes is None:
            otherseqs_bytes = self.get_seqs_bytes_in_pool(pool)
        if do_pool_interact(
            primerpair.all_seq_bytes(),
            otherseqs_bytes,
            self.config.dimer_score,
        ):
            return PrimerPairCheck.INTERACTING

        # Check for other PCR Products
        if detect_new_products(
            primerpair.find_matches(
                self._matchDB,
                remove_expected=False,
                kmersize=self.config.mismatch_kmersize,
                fuzzy=self.config.mismatch_fuzzy,
            ),
            self._matches[pool],
            self.config.mismatch_product_size,
        ):
            return PrimerPairCheck.MISPRIMING

        # If all checks pass, return True
        return PrimerPairCheck.OK

    def setup_coverage(self) -> None:
        """
        Sets up the coverage dict
        :param msa_dict: dict[int, MSA]
        :return: None
        """
        self._coverage = {}
        for msa_index, msa in self._msa_dict.items():
            if msa._mapping_array is None:
                n = len(msa.array[1])
            else:
                n = len(msa._mapping_array)
            self._coverage[msa_index] = np.array([False] * n)

    def recalculate_coverage(self) -> None:
        """
        Recalculates the primertrimmed coverage for all MSA indexes
        :param msa_dict: dict[int, MSA]
        :return: None
        """
        # Create an empty coverage dict
        self.setup_coverage()
        for pp in self.all_primerpairs():
            self.update_coverage(pp, add=True)

    def get_coverage_percent(self, msa_index: int) -> float:
        """
        Returns the coverage percentage for the specified MSA index
        :param msa_index: int
        :return: float | None
        """
        return round(
            self._coverage[msa_index].sum() / len(self._coverage[msa_index]) * 100, 2
        )

    def get_coverage_gaps(self, msa_index: int) -> list[tuple[int, int]]:
        """
        Returns the coverage gaps for the specified MSA index
        :param msa_index: int
        :return: list[tuple[int, int]]
        """
        gaps = []
        in_gap = False
        gap_start = 0
        for i, covered in enumerate(self._coverage[msa_index]):
            if not covered:
                if not in_gap:
                    gap_start = i
                    in_gap = True
            else:
                if in_gap:
                    gaps.append((gap_start, i))
                    in_gap = False
        if in_gap:
            gaps.append((gap_start, i))  # type: ignore
        return gaps

    def update_coverage(self, primerpair: PrimerPair, add: bool = True) -> None:
        """
        Updates the coverage for the specified MSA index
        :param msa_index: int
        :param primerpair: PrimerPair object
        :param add: bool. If True, add to the coverage. If False, remove from the coverage
        :return: None
        """
        # If the msa_index is not in the lookup. Then return
        if primerpair.msa_index not in self._coverage:
            return
        # Check circular
        if primerpair.fprimer.end > primerpair.rprimer.start:
            # Handle circular
            self._coverage[primerpair.msa_index][primerpair.fprimer.region()[1] :] = add
            self._coverage[primerpair.msa_index][: primerpair.rprimer.region()[0]] = (
                add  # slice non inclusive
            )
        else:
            self._coverage[primerpair.msa_index][
                primerpair.fprimer.end : primerpair.rprimer.start
            ] = add

    def next_pool(self) -> int:
        """
        Returns the next pool number.
        Does not directly change self._current_pool
        :return: int
        """
        return (self._current_pool + 1) % self.n_pools

    def add_primer_pair_to_pool(
        self, primerpair: PrimerPair | BedPrimerPair, pool: int, msa_index: int
    ):
        """
        Main method to add a primerpair to a pool. Performs no checks.
        - Adds PrimerPair to the specified pool
        - Updates the PrimerPair's pool and amplicon_number
        - Updates the pools matches
        - Appends PrimerPair to _last_pp_added
        - Sets the Multiplex to the specified pool. Then moves the Multiplex to the next pool


        :param primerpair: PrimerPair object
        :param pool: int
        :param msa_index: int
        :return: None
        """
        # Set the primerpair values
        primerpair.pool = pool
        primerpair.amplicon_number = (
            len(
                [
                    pp
                    for sublist in self._pools
                    for pp in sublist
                    if pp.msa_index == primerpair.msa_index
                ]
            )
            + 1
        )

        # Adds the primerpair's matches to the pools matches
        self._matches[pool].update(
            primerpair.find_matches(
                self._matchDB,
                fuzzy=self.config.mismatch_fuzzy,
                remove_expected=True,
                kmersize=self.config.mismatch_kmersize,
            )
        )

        # Adds the primerpair to the pool
        self._pools[pool].append(primerpair)
        self._current_pool = pool
        self._current_pool = self.next_pool()
        self._last_pp_added.append(primerpair)

        # Update the lookup
        self.update_lookup(primerpair, add=True)
        # Update the coverage
        self.update_coverage(primerpair, add=True)

    def remove_primerpair(self, pp):
        """
        Main method to remove a primerpair from the multiplex
        - Removes the primerpair from the pool
        - Removes the primerpair's matches from the pool's matches
        - Updates the lookup
        - Updates the coverage
        :param pp: PrimerPair object
        :return: None
        """
        # Removes the pp from stack
        self._last_pp_added.remove(pp)

        # Remove the primerpair from the pool
        self._pools[pp.pool].remove(pp)

        # Remove the primerpair's matches from the pool's matches
        self._matches[pp.pool].difference_update(
            pp.find_matches(
                self._matchDB,
                fuzzy=self.config.mismatch_fuzzy,
                remove_expected=False,
                kmersize=self.config.mismatch_kmersize,
            )
        )

        # Update the lookup
        self.update_lookup(pp, add=False)
        # Update the coverage
        self.update_coverage(pp, add=False)

    def remove_last_primer_pair(self) -> PrimerPair:
        """
        This removes the last primerpair added
        - Finds the last primerpair added from self._last_pp_added
        - Removes the primerpair from the pool
        - Removes the primerpair's matches from the pool's matches
        - Moves the current pool to the last primerpair's pool
        - Returns the last primerpair added
        :raises: IndexError if no primerpairs have been added

        :return: PrimerPair object
        """
        last_pp = self._last_pp_added[-1]
        self.remove_primerpair(last_pp)
        self._current_pool = last_pp.pool

        return last_pp

    def does_overlap(self, primerpair: PrimerPair | BedPrimerPair, pool: int) -> bool:
        """
        Does this primerpair overlap with any primerpairs in the pool?
        :param primerpair: PrimerPair object
        :param pool: int
        :return: bool. True if overlaps
        """
        # Circular
        if primerpair.fprimer.end > primerpair.rprimer.start:
            return (
                np.count_nonzero(
                    self._lookup[primerpair.msa_index][
                        pool, primerpair.fprimer.region()[0] :
                    ]
                )
                > 0
                or np.count_nonzero(
                    self._lookup[primerpair.msa_index][
                        pool, : primerpair.rprimer.region()[1]
                    ]
                )
                > 0
            )

        # Get the slice of the lookup
        return (
            np.count_nonzero(
                self._lookup[primerpair.msa_index][
                    pool,
                    primerpair.fprimer.region()[0] : primerpair.rprimer.region()[
                        1
                    ],  # slice non inclusive
                ]
            )
            > 0
        )

    def all_primerpairs(self) -> list[PrimerPair]:
        """
        Returns a list of all primerpairs in the multiplex.
        Sorted by MSA index and amplicon number
        :return: list[PrimerPair]
        """
        all_pp = [pp for pool in (x for x in self._pools) for pp in pool]
        all_pp.sort(key=lambda pp: (str(pp.msa_index), pp.amplicon_number))
        return all_pp

    # Output methods
    def to_bed(
        self,
        headers: list[str] | None = None,
    ) -> str:
        """
        Returns the multiplex as a bed file
        :return: str
        """
        if headers is None:
            headers = ["# artic-bed-version v3.0"]

        return create_bedfile_str(headers, self.all_primerpairs())

    def to_amplicons(
        self,
        trim_primers: bool,
    ) -> str:
        """
        Returns the multiplex as an amplicon file
        :param trim_primers: bool. If True, the primers are trimmed from the amplicons
        :return: str
        """
        return create_amplicon_str(self.all_primerpairs(), trim_primers)

    def polish(
        self,
        msas_dict: dict[int, MSA],
    ) -> None:
        """
        Stochastic optimization to improve the multiplex
        """

        # Create interaction network

        # Find primerpairs that cover uncoverered regions
        msaindex_to_primerpairs = {}
        for msa_index, msa in msas_dict.items():
            for pp in msa.primerpairs:
                numuncoveredpos = (
                    pp.rprimer.start
                    - pp.fprimer.end
                    - np.sum(
                        self._coverage[msa_index][pp.fprimer.end : pp.rprimer.start],
                        dtype=int,
                    )
                )
                if numuncoveredpos > 0:
                    if msa_index not in msaindex_to_primerpairs:
                        msaindex_to_primerpairs[msa_index] = []
                    msaindex_to_primerpairs[msa_index].append((pp, numuncoveredpos))

        pass
