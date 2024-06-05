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
    """

    _pools: list[list[PrimerPair | BedPrimerPair]]
    _current_pool: int
    _last_pp_added: list[PrimerPair]  # Stack to keep track of the last primer added
    _matchDB: MatchDB
    _matches: list[set[tuple]]
    _coverage: dict[int, np.ndarray] | None
    cfg: dict

    def __init__(self, cfg, matchDB: MatchDB) -> None:
        self.n_pools = cfg["npools"]
        self._pools = [[] for _ in range(self.n_pools)]
        self._matches: list[set[tuple]] = [set() for _ in range(self.n_pools)]
        self._current_pool = 0
        self._pp_number = 1
        self.cfg = cfg
        self._matchDB = matchDB
        self._last_pp_added = []
        self._coverage = None

    def setup_coverage(self, msa_dict: dict[int, MSA]) -> None:
        """
        Sets up the coverage dict
        :param n: int. The number of amplicons
        :return: None
        """
        self._coverage = {}
        for msa_index, msa in msa_dict.items():
            if msa._mapping_array is None:
                n = len(msa.array[1])
            else:
                n = len(msa._mapping_array)
            self._coverage[msa_index] = np.array([False] * n)

    def get_coverage_percent(self, msa_index: int) -> float | None:
        if self._coverage is None or msa_index not in self._coverage:
            return None
        return round(
            self._coverage[msa_index].sum() / len(self._coverage[msa_index]) * 100, 2
        )

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

        # Update the coverage
        if self._coverage is not None and msa_index in self._coverage:
            # Check not circular
            if primerpair.start < primerpair.end:
                self._coverage[msa_index][
                    primerpair.fprimer.end : primerpair.rprimer.start
                ] = True
            else:
                # Handle circular
                self._coverage[msa_index][primerpair.fprimer.end :] = True
                self._coverage[msa_index][: primerpair.rprimer.start] = True

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
        # Removes the pp from self._last_pp_added
        last_pp = self._last_pp_added.pop()

        # Remove the primerpair from the pool
        self._pools[last_pp.pool].pop()
        # Remove the primerpair's matches from the pool's matches
        self._matches[last_pp.pool].difference_update(
            last_pp.find_matches(
                self._matchDB,
                fuzzy=self.cfg["mismatch_fuzzy"],
                remove_expected=False,
                kmersize=self.cfg["mismatch_kmersize"],
            )
        )
        # Move the current pool to the last primerpair's pool
        self._current_pool = last_pp.pool

        return last_pp

    def does_overlap(self, primerpair: PrimerPair | BedPrimerPair, pool: int) -> bool:
        """
        Does this primerpair overlap with any primerpairs in the pool?
        :param primerpair: PrimerPair object
        :param pool: int
        :return: bool. True if overlaps
        """
        primerpairs_in_pool = self._pools[pool]

        # Check if the provided primerpair overlaps with any primerpairs in the pool
        for current_primerpairs in primerpairs_in_pool:
            # If they are from the same MSA
            if current_primerpairs.msa_index != primerpair.msa_index:
                # Guard for different MSAs
                continue

            if range(
                max(primerpair.start, current_primerpairs.start),
                min(primerpair.end, current_primerpairs.end) + 1,
            ):
                return True
        # If no overlap
        return False

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
        headers: list[str] | None,
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
