from enum import Enum
from multiprocessing import Pool

# Module imports
from primalscheme3.core.classes import PrimerPair
from primalscheme3.core.multiplex import Multiplex
from primalscheme3.core.mismatches import MatchDB, detect_new_products
from primalscheme3.core.get_window import get_pp_window
from primalscheme3.core.bedfiles import BedPrimerPair
from primalscheme3.core.digestion import _mp_pp_inter_free


# Scheme imports
from primalscheme3.scheme.primer_pair_score import (
    ol_pp_score,
    walk_pp_score,
    bt_ol_pp_score,
)

from primaldimer_py import do_pools_interact_py  # type: ignore


class SchemeReturn(Enum):
    # Added return values
    ADDED_OL_PRIMERPAIR = 1
    ADDED_WALK_PRIMERPAIR = 2
    ADDED_FIRST_PRIMERPAIR = 3
    # Failed return values
    NO_OL_PRIMERPAIR = 4
    NO_WALK_PRIMERPAIR = 5
    NO_FIRST_PRIMERPAIR = 6
    # Misc return values
    ADDED_BACKTRACKED = 7
    NO_BACKTRACK = 8

    ADDED_CIRULAR = 9
    NO_CIRCULAR = 10


class Scheme(Multiplex):
    _pools: list[list[PrimerPair | BedPrimerPair]]
    _current_pool: int
    _last_pp_added: list[PrimerPair]  # Stack to keep track of the last primer added
    _matchDB: MatchDB
    cfg: dict

    def __init__(self, cfg, matchDB: MatchDB):
        self.n_pools = cfg["npools"]
        self._pools = [[] for _ in range(self.n_pools)]
        self._matches: list[set[tuple]] = [set() for _ in range(self.n_pools)]
        self._current_pool = 0
        self._pp_number = 1
        self.cfg = cfg
        self._matchDB = matchDB
        self._last_pp_added = []

    @property
    def npools(self) -> int:
        return self.n_pools

    def add_first_primer_pair(
        self, primerpairs: list[PrimerPair], msa_index
    ) -> SchemeReturn:
        "Adds primerpair to the current pool, and updates the current pool"
        # If there are no primerpairs, return false
        if not primerpairs:
            return SchemeReturn.NO_FIRST_PRIMERPAIR

        # Try and add the first primerpair to an empty pool
        for pool_index in range(self.n_pools):
            if not self._pools[pool_index]:
                self.add_primer_pair_to_pool(primerpairs[0], pool_index, msa_index)
                return SchemeReturn.ADDED_FIRST_PRIMERPAIR

        # Create a hashmap of what seqs are in each pool for quicklook up
        pool_seqs_map: dict[int, list[str]] = {
            index: [
                y
                for sublist in (x.all_seqs() for x in self._pools[index])
                for y in sublist
            ]
            for index in range(self.n_pools)
        }

        # Adds the first valid primerpair
        for primerpair in primerpairs:
            for pool_index in range(self.n_pools):
                if do_pools_interact_py(
                    list(primerpair.all_seqs()),
                    pool_seqs_map[pool_index],
                    self.cfg["dimerscore"],
                ):
                    continue
                if detect_new_products(
                    primerpair.find_matches(
                        self._matchDB,
                        remove_expected=False,
                        kmersize=self.cfg["mismatch_kmersize"],
                        fuzzy=self.cfg["mismatch_fuzzy"],
                    ),
                    self._matches[pool_index],
                    self.cfg["mismatch_product_size"],
                ):
                    continue

                self.add_primer_pair_to_pool(primerpair, pool_index, msa_index)
                return SchemeReturn.ADDED_FIRST_PRIMERPAIR

        # If not primerpair can be added return false
        return SchemeReturn.NO_FIRST_PRIMERPAIR

    def get_leading_coverage_edge(self) -> int:
        """This will return the furthest primer-trimmed region with coverage"""
        # This will crash if no primer has been added, but should not be called until one is
        return self._last_pp_added[-1].rprimer.start

    def get_leading_amplicon_edge(self) -> int:
        """This will return the furthest point of an amplicon"""
        # This will crash if no primer has been added, but should not be called until one is
        return max(self._last_pp_added[-1].rprimer.ends())

    def find_ol_primerpairs(self, all_pp_list, minoverlap) -> list[PrimerPair]:
        """
        Finds all primerpairs that could overlap with the last primerpair added.
        However, it does not check for clash between new PP and the last PP in the same pool
        """
        last_primer_pair: PrimerPair = self._last_pp_added[-1]
        return get_pp_window(
            all_pp_list,
            fp_end_min=last_primer_pair.fprimer.end,
            fp_end_max=last_primer_pair.rprimer.start - minoverlap,
            rp_start_min=max(last_primer_pair.rprimer.ends()) + minoverlap,
        )

    def try_ol_primerpairs(self, all_pp_list, msa_index) -> SchemeReturn:
        """
        This will try and add this primerpair into any valid pool.
        Will return true if the primerpair has been added
        """
        last_pool = self._last_pp_added[-1].pool
        # Find what other pools to look in
        pos_pools_indexes = [
            (last_pool + i) % self.n_pools
            for i in range(self.n_pools)
            if (last_pool + i) % self.n_pools != last_pool
        ]

        # Create a hashmap of all sequences in each pool for quick look up
        index_to_seqs: dict[int, list[str]] = {
            index: [
                y
                for sublist in (x.all_seqs() for x in self._pools[index])
                for y in sublist
            ]
            for index in pos_pools_indexes
        }

        # Find pp that could ol, depending on which pool
        pos_ol_pp = self.find_ol_primerpairs(all_pp_list, self.cfg["min_overlap"])

        # Sort the primerpairs depending on overlap score
        pos_ol_pp.sort(
            key=lambda pp: ol_pp_score(
                pp.rprimer.start,
                len(pp.all_seqs()),
                self.get_leading_coverage_edge() - self.cfg["min_overlap"],
                self.cfg,
            ),
            reverse=True,
        )

        # For each primerpair
        for ol_pp in pos_ol_pp:
            # For each pool
            for pool_index in pos_pools_indexes:
                # If the pool is empty
                if not self._pools[pool_index]:
                    self.add_primer_pair_to_pool(ol_pp, pool_index, msa_index)
                    return SchemeReturn.ADDED_OL_PRIMERPAIR

                # Guard for clash between the last primer in the same pool
                if self._pools[pool_index][-1].msa_index == msa_index and max(
                    self._pools[pool_index][-1].rprimer.ends()
                ) >= min(ol_pp.fprimer.starts()):
                    continue

                # Guard for Primer-Primer Interactions
                if do_pools_interact_py(
                    ol_pp.all_seqs(),
                    index_to_seqs.get(pool_index),
                    self.cfg["dimerscore"],
                ):
                    continue
                # Guard for Primer-Mispriming Products
                if detect_new_products(
                    ol_pp.find_matches(
                        self._matchDB,
                        remove_expected=False,
                        kmersize=self.cfg["mismatch_kmersize"],
                        fuzzy=self.cfg["mismatch_fuzzy"],
                    ),
                    self._matches[pool_index],
                    self.cfg["mismatch_product_size"],
                ):
                    continue

                # If the primer passes all the checks, add it to the pool
                self.add_primer_pair_to_pool(ol_pp, pool_index, msa_index)
                return SchemeReturn.ADDED_OL_PRIMERPAIR

        # If non of the primers work, return false
        return SchemeReturn.NO_OL_PRIMERPAIR

    # backtracking
    def try_backtrack(self, all_pp_list, msa_index) -> SchemeReturn:
        """If there are no other valid ol primerpairs, replace the last primerpair added and try again"""

        # Remove the last primerpair added
        last_pp = self.remove_last_primer_pair()

        # Find all primerpairs that could replace the last primerpair
        pos_ol_pp = [
            pp
            for pp in self.find_ol_primerpairs(all_pp_list, 1)
            if pp != last_pp  # Change minoverlap to 1 to help the solver
        ]
        # Sort the primerpair on score
        pos_ol_pp.sort(
            key=lambda pp: bt_ol_pp_score(
                pp.rprimer.start,
                len(pp.all_seqs()),
                self.get_leading_coverage_edge() - 1,
                self.cfg,
            ),
            reverse=True,
        )

        # Find what other pools to look in
        pos_pools_indexes = [last_pp.pool]
        # Create a hashmap of all sequences in each pool for quick look up
        index_to_seqs: dict[int, list[str]] = {
            index: [
                y
                for sublist in (x.all_seqs() for x in self._pools[index])
                for y in sublist
            ]
            for index in pos_pools_indexes
        }

        # For each replacement primerpair
        for ol_pp in pos_ol_pp:
            # For each pool
            for pool_index in pos_pools_indexes:
                # If the pool is empty
                if not self._pools[pool_index]:
                    self.add_primer_pair_to_pool(ol_pp, pool_index, msa_index)
                    return SchemeReturn.ADDED_BACKTRACKED

                # If the last primer is from the same msa and does clash, skip it
                if self._pools[pool_index][-1].msa_index == msa_index and max(
                    self._pools[pool_index][-1].rprimer.ends()
                ) >= min(ol_pp.fprimer.starts()):
                    continue

                # Guard for interactions
                if do_pools_interact_py(
                    ol_pp.all_seqs(),
                    index_to_seqs.get(pool_index),
                    self.cfg["dimerscore"],
                ):
                    continue

                # Guard for misprimal
                if detect_new_products(
                    ol_pp.find_matches(
                        self._matchDB,
                        remove_expected=False,
                        kmersize=self.cfg["mismatch_kmersize"],
                        fuzzy=self.cfg["mismatch_fuzzy"],
                    ),
                    self._matches[pool_index],
                    self.cfg["mismatch_product_size"],
                ):
                    continue
                # If all checks pass add the primerpair to the pool this is a valid alternative
                # See if the valid alternative had a valid overlap

                # Add the primer to the pool
                self.add_primer_pair_to_pool(ol_pp, pool_index, msa_index)
                # Try and add an overlap
                ol_result = self.try_ol_primerpairs(all_pp_list, msa_index=msa_index)
                if ol_result == SchemeReturn.ADDED_OL_PRIMERPAIR:
                    # Fixed the problem
                    return SchemeReturn.ADDED_BACKTRACKED
                else:
                    # No overlap was found. Remove the primerpair and the try the next one
                    self.remove_last_primer_pair()
                    continue

        # If non of the primers work, add the last pp back in and return false
        self.add_primer_pair_to_pool(last_pp, last_pp.pool, msa_index)
        return SchemeReturn.NO_BACKTRACK

    def try_walk_primerpair(self, all_pp_list, msa_index) -> SchemeReturn:
        """
        Find the next valid primerpair while walking forwards
        """
        last_pool = self._last_pp_added[-1].pool
        # Find what other pools to look in, can look in same pool
        pos_pools_indexes = [
            (last_pool + i) % self.n_pools for i in range(self.n_pools)
        ]

        # Create a hashmap of all sequences in each pool for quick look up
        index_to_seqs: dict[int, list[str]] = {
            index: [
                y
                for sublist in (x.all_seqs() for x in self._pools[index])
                for y in sublist
            ]
            for index in pos_pools_indexes
        }
        # Find the walking start index
        walking_min = self._last_pp_added[-1].rprimer.start - self.cfg["min_overlap"]

        # Find the first primer that could walk
        ## Use that index to slice the list
        pos_walk_pp = []
        for index, pp in enumerate(all_pp_list):
            if pp.fprimer.end > walking_min:
                pos_walk_pp = all_pp_list[index:]
                break

        # Sort walking primers by score
        pos_walk_pp.sort(
            key=lambda pp: walk_pp_score(
                pp.fprimer.end, len(pp.all_seqs()), self._last_pp_added[-1].end
            ),
            reverse=True,
        )

        # For each primer, try each pool
        for walk_pp in pos_walk_pp:
            for pool_index in pos_pools_indexes:
                # If the pool is empty add the first primer
                if not self._pools[pool_index]:
                    self.add_primer_pair_to_pool(walk_pp, pool_index, msa_index)
                    return SchemeReturn.ADDED_WALK_PRIMERPAIR

                # Guard for clash between the last primer in the same pool
                if self._pools[pool_index][-1].msa_index == msa_index and max(
                    self._pools[pool_index][-1].rprimer.ends()
                ) >= min(walk_pp.fprimer.starts()):
                    continue

                # Guard for Primer-Primer Interactions
                if do_pools_interact_py(
                    walk_pp.all_seqs(),
                    index_to_seqs.get(pool_index),
                    self.cfg["dimerscore"],
                ):
                    continue
                # Guard for Primer-Mispriming Products
                if detect_new_products(
                    walk_pp.find_matches(
                        self._matchDB,
                        remove_expected=False,
                        kmersize=self.cfg["mismatch_kmersize"],
                        fuzzy=self.cfg["mismatch_fuzzy"],
                    ),
                    self._matches[pool_index],
                    self.cfg["mismatch_product_size"],
                ):
                    continue
                # If the primer passes all the checks, add it to the pool
                self.add_primer_pair_to_pool(walk_pp, pool_index, msa_index)
                return SchemeReturn.ADDED_WALK_PRIMERPAIR
        # If non of the primers work, return false
        return SchemeReturn.NO_WALK_PRIMERPAIR

    def try_circular(self, msa):
        """
        This will try and add a primerpair that can span from the end of the msa back to the start as if the genome was circular
        """
        first_pp: PrimerPair = self._last_pp_added[0]
        last_pp: PrimerPair = self._last_pp_added[-1]
        last_pool = last_pp.pool

        # Find all possible fkmers and rkmer that could span the end of the msa
        pos_fkmers = [
            fkmer
            for fkmer in msa.fkmers
            if fkmer.end < last_pp.rprimer.start
            and fkmer.end > last_pp.rprimer.start - 200
        ]
        pos_rkmers = [
            rkmer
            for rkmer in msa.rkmers
            if rkmer.start > first_pp.fprimer.end
            and rkmer.start < first_pp.fprimer.end + 200
        ]

        # Create all the primerpairs
        non_checked_pp = []
        for fkmer in pos_fkmers:
            for rkmer in pos_rkmers:
                non_checked_pp.append(PrimerPair(fkmer, rkmer, msa.msa_index))

        ## Interaction check all the primerpairs
        iter_free_primer_pairs = []
        with Pool(self.cfg["n_cores"]) as p:
            mp_pp_bool = p.imap_unordered(
                _mp_pp_inter_free, ((pp, self.cfg) for pp in non_checked_pp)
            )
            for bool, pp in zip(mp_pp_bool, non_checked_pp):
                if not bool:
                    iter_free_primer_pairs.append(pp)

        # Sort the primerpairs by number of primers
        iter_free_primer_pairs.sort(key=lambda pp: len(pp.all_seqs()))

        pos_pools_indexes = [
            (last_pool + i) % self.n_pools for i in range(self.n_pools)
        ]

        # Create a hashmap of all sequences in each pool for quick look up
        index_to_seqs: dict[int, list[str]] = {
            index: [
                y
                for sublist in (x.all_seqs() for x in self._pools[index])
                for y in sublist
            ]
            for index in pos_pools_indexes
        }

        for c_pp in iter_free_primer_pairs:
            for pool_index in pos_pools_indexes:
                # Guard for clash between the last primer in the same pool
                if self._pools[pool_index][-1].msa_index == msa.msa_index and max(
                    self._pools[pool_index][-1].rprimer.ends()
                ) >= min(c_pp.fprimer.starts()):
                    continue

                # Guard for clash between the first primer in the same pool
                if self._pools[pool_index][0].msa_index == msa.msa_index and min(
                    self._pools[pool_index][0].fprimer.starts()
                ) <= max(c_pp.rprimer.ends()):
                    continue

                # Guard for Primer-Primer Interactions
                if do_pools_interact_py(
                    c_pp.all_seqs(),
                    index_to_seqs.get(pool_index),
                    self.cfg["dimerscore"],
                ):
                    continue

                # Skip the miss priming product check, as we are using special indexes
                # If the primer passes all the checks, add it to the pool
                self.add_primer_pair_to_pool(c_pp, pool_index, msa.msa_index)
                return SchemeReturn.ADDED_CIRULAR

        # No Primers could be added
        return SchemeReturn.NO_CIRCULAR

    def all_primers(self) -> list[PrimerPair]:
        all_pp = [pp for pool in (x for x in self._pools) for pp in pool]
        all_pp.sort(key=lambda pp: (str(pp.msa_index), pp.amplicon_number))
        return all_pp