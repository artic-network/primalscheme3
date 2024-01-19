from primalscheme3.core.classes import PrimerPair, MatchDB
from primalscheme3.core.bedfiles import BedPrimerPair


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
        primerpair.amplicon_number = len(
            [
                pp
                for sublist in self._pools
                for pp in sublist
                if pp.msa_index == primerpair.msa_index
            ]
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

    def remove_last_primer_pair(self) -> PrimerPair:
        """
        This removes the last primerpair added
        - Finds the last primerpair added from self._last_pp_added
        - Removes the primerpair from the pool
        - Removes the primerpair's matches from the pool's matches
        - Moves the current pool to the last primerpair's pool
        - Returns the last primerpair added

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
        headers: list[str] | None = ["# artic-bed-version v3.0"],
    ) -> str:
        """
        Returns the multiplex as a bed file
        :return: str
        """
        primer_bed_str: list[str] = []

        # Ensure headers are commented and valid
        if headers is not None:
            for headerline in headers:
                if not headerline.startswith("#"):
                    headerline = "# " + headerline
                primer_bed_str.append(headerline.strip())

        # Add the primerpairs to the bed file
        for pp in self.all_primerpairs():
            primer_bed_str.append(pp.to_bed().strip())

        return "\n".join(primer_bed_str)
