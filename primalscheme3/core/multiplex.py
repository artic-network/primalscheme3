import numpy as np

from primalscheme3.core.bedfiles import (
    BedPrimerPair,
    create_amplicon_str,
    create_bedfile_str,
)
from primalscheme3.core.classes import MatchDB, PrimerPair
from primalscheme3.core.msa import MSA


class Multiplex:
    """
    This is the baseclass for all multiplexes (Scheme / Panel)
    - It allows mutliple pools

    Params:
    - cfg: dict. The configuration dict
    - matchDB: MatchDB object
    - msa_dict: dict[int, MSA]. A dict of MSA objects

    Internal:
    - _pools: list[list[PrimerPair | BedPrimerPair]]. List of pools with PrimerPairs
    - _current_pool: int. The current pool number
    - _last_pp_added: list[PrimerPair]. Stack to keep track of the last primer added
    - _matchDB: MatchDB object
    - _matches: list[set[tuple]]. A list of sets with matches for each pool
    - _msa_dict: dict[int, MSA]. A dict of MSA objects
    - _coverage: dict[int, np.ndarray]. PrimerTrimmed Coverage for each MSA index
    - _lookup: dict[int, np.ndarray] | None. A lookup for the primerpairs in the multiplex
    """

    _pools: list[list[PrimerPair | BedPrimerPair]]
    _current_pool: int
    _last_pp_added: list[PrimerPair]  # Stack to keep track of the last primer added
    _matchDB: MatchDB
    _matches: list[set[tuple]]
    _msa_dict: dict[int, MSA]
    _coverage: dict[int, np.ndarray]  # PrimerTrimmed Coverage for each MSA index
    _lookup: dict[int, np.ndarray]
    cfg: dict

    def __init__(self, cfg, matchDB: MatchDB, msa_dict: dict[int, MSA]) -> None:
        self.n_pools = cfg["npools"]
        self._pools = [[] for _ in range(self.n_pools)]
        self._matches: list[set[tuple]] = [set() for _ in range(self.n_pools)]
        self._current_pool = 0
        self._pp_number = 1
        self.cfg = cfg
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
        circular = primerpair.start >= primerpair.end

        # If the msa_index is not in the lookup. Then return
        if primerpair.msa_index not in self._lookup:
            return

        if add:
            value = primerpair
        else:
            value = None

        if circular:
            self._lookup[primerpair.msa_index][primerpair.pool, primerpair.start :] = (
                value
            )
            self._lookup[primerpair.msa_index][primerpair.pool, : primerpair.end] = (
                value
            )
        else:
            self._lookup[primerpair.msa_index][
                primerpair.pool, primerpair.start : primerpair.end
            ] = value

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

    def calculate_coverage(self) -> None:
        """
        Recalculates the coverage for all MSA indexes
        :param msa_dict: dict[int, MSA]
        :return: None
        """
        for pp in self.all_primerpairs():
            # Skip PrimerPairs with no MSA index
            if pp.msa_index not in self._coverage:
                continue

            # Check not circular
            if pp.start < pp.end:
                self._coverage[pp.msa_index][pp.fprimer.end : pp.rprimer.start] = True
            else:
                # Handle circular
                self._coverage[pp.msa_index][pp.fprimer.end :] = True
                self._coverage[pp.msa_index][: pp.rprimer.start] = True

    def get_coverage_percent(self, msa_index: int) -> float:
        """
        Returns the coverage percentage for the spesified MSA index
        :param msa_index: int
        :return: float | None
        """
        return round(
            self._coverage[msa_index].sum() / len(self._coverage[msa_index]) * 100, 2
        )

    def update_coverage(
        self, msa_index: int, primerpair: PrimerPair, add: bool = True
    ) -> None:
        """
        Updates the coverage for the spesified MSA index
        :param msa_index: int
        :param primerpair: PrimerPair object
        :param add: bool. If True, add to the coverage. If False, remove from the coverage
        :return: None
        """
        # If the msa_index is not in the lookup. Then return
        if primerpair.msa_index not in self._coverage:
            return
        # Check not circular
        if primerpair.start < primerpair.end:
            self._coverage[msa_index][
                primerpair.fprimer.end : primerpair.rprimer.start
            ] = add
        else:
            # Handle circular
            self._coverage[msa_index][primerpair.fprimer.end :] = add
            self._coverage[msa_index][: primerpair.rprimer.start] = add

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
        - Adds PrimerPair to the spesified pool
        - Updates the PrimerPair's pool and amplicon_number
        - Updates the pools matches
        - Appends PrimerPair to _last_pp_added
        - Sets the Mutliplex to the spesified pool. Then moves the Mutliplex to the next pool


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
                fuzzy=self.cfg["mismatch_fuzzy"],
                remove_expected=True,
                kmersize=self.cfg["mismatch_kmersize"],
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
        self.update_coverage(msa_index, primerpair, add=True)

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
                fuzzy=self.cfg["mismatch_fuzzy"],
                remove_expected=False,
                kmersize=self.cfg["mismatch_kmersize"],
            )
        )

        # Update the lookup
        self.update_lookup(pp, add=False)
        # Update the coverage
        self.update_coverage(pp.msa_index, pp, add=False)

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
        # Get the slice of the lookup
        lookup_slice = self._lookup[primerpair.msa_index][pool][
            primerpair.start : primerpair.end - 1
        ]
        return np.count_nonzero(lookup_slice) > 0

    def all_primerpairs(self) -> list[PrimerPair]:
        """
        Returns a list of all primerpairs in the multiplex.
        Sorted by MSA index and amplicon number
        :return: list[PrimerPair]
        """
        all_pp = [pp for pool in (x for x in self._pools) for pp in pool]
        all_pp.sort(key=lambda pp: (str(pp.msa_index), pp.amplicon_number))
        return all_pp

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
