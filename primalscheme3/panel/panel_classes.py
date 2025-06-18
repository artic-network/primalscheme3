import pathlib
from enum import Enum
from itertools import islice

import numpy as np

# Iterations checker
from primalscheme3.core.bedfiles import BedPrimerPair

# Core Module imports
from primalscheme3.core.classes import PrimerPair
from primalscheme3.core.config import Config
from primalscheme3.core.errors import CustomErrors
from primalscheme3.core.mismatches import MatchDB
from primalscheme3.core.msa import MSA
from primalscheme3.core.multiplex import Multiplex, PrimerPairCheck
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.core.seq_functions import (
    entropy_score_array,
)


class PanelRunModes(Enum):
    ENTROPY = "entropy"
    REGION_ONLY = "region-only"
    EQUAL = "equal"


class PanelReturn(Enum):
    """
    Enum for the return values of the Panel class.
    """

    ADDED_PRIMERPAIR = 0  # A primerpair was added
    NO_PRIMERPAIRS = 1  # No more primerpairs to add (finished)
    NO_PRIMERPAIRS_IN_MSA = 3  # No more primerpairs in the current MSA (moving to next)


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
    group: str | None

    def __init__(
        self,
        chromname: str,
        start: int,
        stop: int,
        name: str,
        score: int,
        group: str | None = None,
    ) -> None:
        self.chromname = str(chromname).strip()
        self.start = int(start)
        self.stop = int(stop)
        self.name = str(name).strip()
        self.score = int(score)
        if self.start >= self.stop:
            raise ValueError(f"{self.name}: Circular regions are not supported.")

        # Parse Group
        if group is not None and group != "":
            group = str(group).strip()
        else:
            group = None

        self.group = group

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
        group_str = self.group if self.group else ""
        return f"{self.chromname}\t{self.start}\t{self.stop}\t{self.name}\t{self.score}\t{group_str}\n"


class RegionParser:
    @staticmethod
    def from_list(bed_list: list[str]) -> Region:
        """
        Returns a single region object from a list representing a line of a bedfile, split by tabs
        """
        # Handle extra columns
        if len(bed_list) >= 5:
            chromname, start, stop, name, score = bed_list[:5]
        else:
            raise CustomErrors(
                f"Invalid region: {' '.join(bed_list)}. Requires 5 or more columns."
            )

        try:
            start = int(start)
            stop = int(stop)
            score = int(score)
        except ValueError as e:
            raise CustomErrors(
                f"Invalid region: {' '.join(bed_list)}. Start, stop, score must be integers."
            ) from e

        # Handle group column
        try:
            group = bed_list[5]
            if group == "" or not group:  # Allow empty group
                group = None
        except IndexError:
            group = None

        return Region(chromname, start, stop, name, score, group)

    @staticmethod
    def from_str(bed_str: str) -> Region:
        """
        Returns a single region object from a string representing a line of a bedfile
        """
        return RegionParser.from_list(bed_str.strip().split("\t"))


class PanelMSA(MSA):
    # Score arrays
    _score_array: np.ndarray | None
    _midx_entropy_array: np.ndarray | None
    regions: list[Region] | None
    region_group_count: dict[str, int] | None
    primerpairpointer: int

    def __init__(
        self,
        name: str,
        path: pathlib.Path,
        msa_index: int,
        mapping: str,
        progress_manager: ProgressManager,
        config: Config,
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
            config=config,
        )

        # Create the primerpairpointer
        self.primerpairpointer = 0
        self.regions = None

    def create_score_array(
        self,
        regions: list[Region] | None,
        mode: PanelRunModes = PanelRunModes.REGION_ONLY,
    ) -> None:
        match mode:
            case PanelRunModes.REGION_ONLY:
                if regions is None:
                    raise ValueError("Regions must be provided to create score array.")
                self._score_array = np.zeros(len(self._mapping_array), dtype=int)
                for region in regions:
                    self._score_array[region.start : region.stop] += region.score
            case PanelRunModes.ENTROPY:
                # Turn the msa indexed entropy array into a ref indexed array
                msa_entropy_array = entropy_score_array(self.array)
                self._score_array = np.array(
                    [
                        x + 0.01  # Add a small value to avoid 0
                        for msai, x in enumerate(msa_entropy_array)
                        if self._mapping_array[msai] is not None
                    ]
                )
            case PanelRunModes.EQUAL:
                self._score_array = np.ones(len(self._mapping_array), dtype=int)

    def create_entropy_array(self):
        self._midx_entropy_array = np.array(entropy_score_array(self.array))

    def add_regions(self, regions: list[Region] | None) -> None:
        self.regions = regions
        if regions is not None:
            # Create an empty dict
            region_groups = {x.group for x in regions if x.group}
            self.region_group_count = {x: 0 for x in region_groups}

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
        assert self._midx_entropy_array is not None
        ppptstart, ppptend = pp.primertrimmed_region()

        return np.sum(
            self._midx_entropy_array[
                self._ref_to_msa[ppptstart] : self._ref_to_msa[ppptend]
            ]
        )

    def get_pp_score(self, pp: PrimerPair) -> int:
        """
        Returns number of SNPs in the primertrimmed amplicon
        """
        assert self._score_array is not None
        ppptstart, ppptend = pp.primertrimmed_region()

        # circular
        if pp.fprimer.end > pp.rprimer.start:
            return np.sum(self._score_array[ppptstart:]) + np.sum(
                self._score_array[:ppptend]
            )
        else:
            return np.sum(self._score_array[ppptstart:ppptend])

    def update_score_array(self, addedpp: PrimerPair, newscore: int = 0) -> None:
        """
        Updates the score array with the added primerpair
        """
        assert self._score_array is not None
        ppptstart, ppptend = addedpp.primertrimmed_region()
        # circular
        if addedpp.fprimer.end > addedpp.rprimer.start:
            self._score_array[ppptstart:] = newscore
            self._score_array[:ppptend] = newscore
        else:
            self._score_array[ppptstart:ppptend] = newscore


class Panel(Multiplex):
    # Base class
    _pools: list[list[PrimerPair | BedPrimerPair]]
    _current_pool: int
    _last_pp_added: list[PrimerPair]  # Stack to keep track of the last primer added
    _matchDB: MatchDB
    config: Config
    _msa_dict: dict[int, PanelMSA]  # type: ignore

    # New attributes
    _current_msa_index: int

    # Keep adding
    _is_msa_index_finished: dict[int, bool]

    def __init__(
        self,
        msa_dict: dict[int, PanelMSA],
        config: Config,
        matchdb: MatchDB,
        logger=None,
    ) -> None:
        super().__init__(config, matchdb, msa_dict)  # type: ignore

        self._current_msa_index = 0
        self._failed_primerpairs = [set() for _ in range(self.n_pools)]

        # Keep adding
        self._is_msa_index_finished = {msa_index: False for msa_index in msa_dict}

        self.logger = logger

    def _next_msa(self) -> int:
        """
        Updates the current msa index to the next msa. Returns the new msa index.
        :return: The new msa index.
        """
        self._current_msa_index = (self._current_msa_index + 1) % len(self._msa_dict)
        return self._current_msa_index

    def _add_primerpair(
        self, primerpair: PrimerPair, pool: int, msa_index: int
    ) -> None:
        # Add the primerpair to the pool
        super().add_primer_pair_to_pool(primerpair, pool, msa_index)

        # Update the msa Score array
        if msa_index in self._msa_dict:
            self._msa_dict[msa_index].update_score_array(primerpair)

        # Update the current MSA
        self._next_msa()

    def add_next_primerpair(
        self, max_amplicons_group: None | int = None
    ) -> PanelReturn:
        """
        Try and add the next primerpair
        """
        # Get current MSA
        current_msa = self._msa_dict[self._current_msa_index]
        current_pool = self._current_pool

        # Check if the current MSA is finished
        if self._is_msa_index_finished[self._current_msa_index]:
            self._next_msa()
            return PanelReturn.NO_PRIMERPAIRS_IN_MSA

        # Pools to try
        pos_pools_indexes = [
            (current_pool + i) % self.n_pools for i in range(self.n_pools)
        ]

        # All seqs in each pool
        seqs_in_each_pool = {
            pos_pool: [
                seq
                for seq in (pp.all_seq_bytes() for pp in self._pools[pos_pool])
                for seq in seq
            ]
            for pos_pool in pos_pools_indexes
        }

        # Remove primerpairs with no score
        # current_msa.primerpairs = [
        #     pp for pp in current_msa.primerpairs if current_msa.get_pp_score(pp) > 0
        # ]
        # For the primerpairs in the current MSA Sort all primerpairs by score
        current_msa.primerpairs.sort(
            key=lambda x: (
                self.calc_pp_score(current_msa, x),
                -len(x.all_seqs()),
            ),
            reverse=True,
        )
        for pospp in current_msa.primerpairs:
            # Check primer has score
            if self.calc_pp_score(current_msa, pospp, False) <= 0:
                continue
            for pospool in pos_pools_indexes:
                # Check if the primerpair can be added
                match self.check_primerpair_can_be_added(
                    pospp, pospool, seqs_in_each_pool[pospool]
                ):
                    case PrimerPairCheck.OK:
                        # _add_primerpair updates the score array for the added primerpair
                        self._add_primerpair(pospp, pospool, self._current_msa_index)
                        # Check for other regions in the same group
                        # Find regions with overlap
                        if current_msa.regions:
                            # Find other regions in the same group
                            other_regions = self.get_regions_with_group_overlap()
                            # Update the MSA region_group_count
                            for region in other_regions:
                                assert current_msa.region_group_count is not None
                                assert current_msa._score_array is not None
                                assert region.group is not None

                                if max_amplicons_group is not None and (
                                    current_msa.region_group_count[region.group]
                                    >= max_amplicons_group
                                ):  # type: ignore
                                    current_msa._score_array[
                                        region.start : region.stop
                                    ] = 0

                                current_msa.region_group_count[region.group] += 1  # type: ignore

                        return PanelReturn.ADDED_PRIMERPAIR
                    case _:
                        continue
        # No more primerpairs in the current MSA
        self._is_msa_index_finished[self._current_msa_index] = True
        return PanelReturn.NO_PRIMERPAIRS_IN_MSA

    def calc_pp_score(self, msa: PanelMSA, pp: PrimerPair, bonus_ol=True) -> int:
        """
        Returns number of SNPs in the primertrimmed amplicon
        """
        # Get the score
        score = msa.get_pp_score(pp)

        # Apply a bonus for good gc content
        gc_diff = pp.get_score()  # 0-0.5
        score = score * (1 - gc_diff)

        # Apply a bonus if the primerpair spans a coverage break
        if bonus_ol:
            coverage_slice = self._coverage[msa.msa_index][
                pp.fprimer.end : pp.rprimer.start
            ]
            if True in coverage_slice and False in coverage_slice:
                score = score * 10

        return int(score)

    def get_regions_with_group_overlap(self) -> list[Region]:
        """
        Finds all regions groups which are overlapping with the last primerpair added.
        Returns all regions with the same groups.
        """
        last_primerpair = self._last_pp_added[-1]
        current_msa = self._msa_dict[self._current_msa_index]
        assert current_msa.regions is not None

        ol_regions_groups = {
            x.group
            for x in current_msa.regions
            if x.group
            and does_overlap(
                (
                    last_primerpair.fprimer.region()[0],
                    last_primerpair.rprimer.region()[1],
                    self._current_msa_index,
                ),
                [(x.start, x.stop, self._current_msa_index)],
            )
        }
        other_regions = [x for x in current_msa.regions if x.group in ol_regions_groups]
        return other_regions
