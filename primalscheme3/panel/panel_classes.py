# Core imports
from primalscheme3.core.classes import FKmer, RKmer, PrimerPair
from primalscheme3.core.digestion import digest
from primalscheme3.core.get_window import get_r_window_FAST2

# from score_function import check_extention_seqs, check_kmers, check_seqs
from primaldimer_py import do_pools_interact_py  # type: ignore

import random
from math import exp, sqrt
from numpy import ndarray
import networkx as nx
import itertools


def should_confirm_change(old_score, new_score, g=1) -> bool:
    """Returns True if the change should be confirmed"""

    # A negative prop_change means a decrease in score
    prop_change = (new_score - old_score) / old_score

    # If improvement or acceptable loss
    if prop_change < 0 or random.random() < exp(-sqrt(prop_change)) / g:
        return True
    else:
        return False


class Region:
    """
    This is the main class that contains region (genome infomation)
    """

    ref_name: str
    start: int
    end: int
    msa: ndarray

    fkmers: list[FKmer]
    rkmers: list[RKmer]

    primer_pairs: list[PrimerPair]
    _single_primer_pairs: list[PrimerPair]

    current_primer_pair: PrimerPair
    _tmp_primer_pair: PrimerPair

    def __init__(self, bedline) -> None:
        self.ref_name = bedline[0].strip()
        self.start = int(bedline[1])
        self.end = int(bedline[2])

        self._single_primer_pairs = None

    def add_msa(self, msa: ndarray, amplicon_max_length: int) -> None:
        "This will extract the region from the msa + one amplicon length to either side"
        self.msa = msa[
            :, self.start - amplicon_max_length : self.end + amplicon_max_length
        ]

    def digest(self, cfg, thermo_cfg):
        "Generate the all the posaible fkmer and rkmers from the msa slice"
        self.fkmers, self.rkmers = digest(
            msa_array=self.msa,
            cfg=cfg,
            thermo_cfg=thermo_cfg,
            offset=self.start
            - cfg[
                "amplicon_size_max"
            ],  # As the msa is a slice, we need to add an offset to ensure primer index is the genomic index rather than the slice index
        )

    def generate_all_primerpairs(self, cfg):
        primer_pairs: list[PrimerPair] = []

        for f in self.fkmers:
            pos_r = get_r_window_FAST2(
                kmers=self.rkmers,
                start=min(f.starts()) + cfg["amplicon_size_min"],
                end=min(f.starts()) + cfg["amplicon_size_max"],
            )
            for r in pos_r:
                primer_pairs.append(PrimerPair(f, r))
        self.primer_pairs = primer_pairs

    def find_single_pp_coverage(self, cfg) -> list[PrimerPair]:
        "Will find all single primerpairs that cover the whole region or return empty list"
        if self._single_primer_pairs is None:
            single_pp = []
            # Check if primerpairs can span the region
            if self.end - self.start > cfg["amplicon_size_max"]:
                pass
            else:
                for pp in self.primer_pairs:
                    if (
                        pp.fprimer.end < self.start
                        and pp.rprimer.start > self.end
                        and not do_pools_interact_py(
                            list(pp.fprimer.seqs),
                            list(pp.rprimer.seqs),
                            cfg["dimer_threshold"],
                        )
                    ):
                        single_pp.append(pp)

            # If there are amplcions asign the best
            if single_pp:
                single_pp.sort(key=lambda pp: len(pp.all_seqs()))
                self.current_primer_pair = single_pp[0]

            self._single_primer_pairs = single_pp
            return single_pp

        else:
            return self._single_primer_pairs

    def __hash__(self) -> int:
        return hash(f"{self.ref_name}*{self.start}*{self.end}")

    def __eq__(self, other):
        if isinstance(other, Region):
            return self.__hash__() == other.__hash__()
        else:
            return False

    def __str__(self) -> str:
        return f"{self.ref_name}-{self.start}-{self.end}"

    def test_new_amplicon(self) -> None:
        """Picks a random amplicon for testing, weighted towards less sequences"""
        self._tmp_primer_pair = random.choices(
            self._single_primer_pairs,
            [1 / len(x.all_seqs()) for x in self._single_primer_pairs],
        )[0]

    def confirm_new_amplicon(self) -> None:
        """Replaced the confirmed amplicons with the tmp amplicon"""
        self.current_primer_pair = self._tmp_primer_pair

    def get_amplicon(self) -> PrimerPair:
        """Returns the confirmed amplicon"""
        if self.current_primer_pair is None:
            self.current_primer_pair = random.choice(self.amplicon_options)
            return self.current_primer_pair
        else:
            return self.current_primer_pair


class Scheme:
    regions: list[Region]
    all_primer_seqs: list[str]

    edge_map: nx.Graph

    def __init__(self, regions: list[Region]) -> None:
        g = nx.Graph()

        # Add all nodes
        g.add_nodes_from(regions)

        # Add all edges between regions if there is an interaction
        for r1, r2 in itertools.product(regions, regions):
            if r1 != r2 and not do_pools_interact_py(
                list(r1.current_primer_pair.all_seqs()),
                list(r2.current_primer_pair.all_seqs()),
                -22,
            ):
                g.add_edge(r1, r2)

        self.edge_map = g

    def remove_primerpair(self, region) -> None:
        """
        This will remove all egdes linked to a region
        """
        self.edge_map.remove_node(region)
        self.edge_map.add_node(region)

    def add_primerpair(self, region) -> None:
        """
        This will add edges betweem the regions current primerpair and all other regions
        """
        for old_region in self.edge_map.nodes:
            if region != old_region and not do_pools_interact_py(
                list(region.current_primer_pair.all_seqs()),
                list(old_region.current_primer_pair.all_seqs()),
                -22,
            ):
                self.edge_map.add_edge(region, old_region)

    def total_interactions(self) -> int:
        "Return the total number of interactions"
        return self.edge_map.number_of_edges()

    def interactions_map(self) -> dict[Region:int]:
        "Returns a dict of how many interactions each region has"
        return {x[0]: len(x[1]) for x in self.edge_map.adjacency()}

    def find_worst_region(self) -> Region:
        "Returns the region with the max number of interactions, if draw might be undefined"
        data = self.interactions_map()
        return max(data, key=data.get)

    def replace_worst_primerpair(self) -> None:
        "This will pick a bad reigon, select a new amplicon from the posiable options, and then recalculate the score"
        region = self.find_worst_region()

        # Removes the old primerpair
        self.remove_primerpair(region)

        # This picks a new amplicon and then asigns it to the confirmed spot
        region.test_new_amplicon()
        region.confirm_new_amplicon()

        # Calculates the new score with the new amplicon
        self.add_primerpair(region)

    def to_bed_format(self) -> str:
        "Returns all primers in the primer.bed format"
        regions: list[Region] = list(self.edge_map.nodes)
        regions.sort(key=lambda r: (r.ref_name, r.start))
        return "".join(
            [r.current_primer_pair.__str__(r.ref_name, "panel") for r in regions]
        )
