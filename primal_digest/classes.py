from primal_digest.config import AMBIGUOUS_DNA_COMPLEMENT
from primal_digest.iteraction import all_inter_checker
from primal_digest.get_window import get_pp_window

class FKmer:
    end: int
    seqs: set[str]

    def __init__(self, end, seqs) -> None:
        self.end = end
        self.seqs = seqs

    def len(self) -> set[int]:
        return {len(x) for x in self.seqs}

    def starts(self) -> set[int]:
        return {self.end - len(x) for x in self.seqs}

    def __str__(self, referance, amplicon_prefix, pool) -> str:
        string_list = []
        counter = 0
        seqs = list(self.seqs)
        seqs.sort()
        for seq in seqs:
            string_list.append(
                f"{referance}\t{self.end-len(seq)}\t{self.end}\t{amplicon_prefix}_LEFT_{counter}\t{pool}\t+\t{seq}\n"
            )
            counter += 1
        return "".join(string_list)

    def __hash__(self) -> int:
        seqs = list(self.seqs)
        seqs.sort()
        return hash(f"{self.end}{self.seqs}")

    def __eq__(self, other):
        if isinstance(other, FKmer):
            return self.__hash__() == other.__hash__()


class RKmer:
    start: int
    seqs: set[str]

    def __init__(self, start, seqs) -> None:
        self.start = start
        self.seqs = seqs

    def len(self) -> set[int]:
        return {len(x) for x in self.seqs}

    def ends(self) -> set[int]:
        return {len(x) + self.start for x in self.seqs}

    def __str__(self, referance, amplicon_prefix, pool) -> str:
        string_list = []
        counter = 0
        seqs = list(self.seqs)
        seqs.sort()
        for seq in seqs:
            string_list.append(
                f"{referance}\t{self.start}\t{self.start+len(seq)}\t{amplicon_prefix}_RIGHT_{counter}\t{pool}\t-\t{seq}\n"
            )
            counter += 1
        return "".join(string_list)

    def reverse_complement(self):
        rev_seqs = {x[::-1] for x in self.seqs}
        self.seqs = {
            "".join(AMBIGUOUS_DNA_COMPLEMENT[base.upper()] for base in seq)
            for seq in rev_seqs
        }
        return self

    def __hash__(self) -> int:
        seqs = list(self.seqs)
        seqs.sort()
        return hash(f"{self.start}{self.seqs}")

    def __eq__(self, other):
        if isinstance(other, RKmer):
            return self.__hash__() == other.__hash__()


class PrimerPair:
    fprimer: FKmer
    rprimer: RKmer
    amplicon_number: int
    pool: int
    msa_index: int

    def __init__(self, fprimer, rprimer, amplicon_number=-1, pool=-1, msa_index=-1):
        self.fprimer = fprimer
        self.rprimer = rprimer
        self.amplicon_number = amplicon_number
        self.pool = pool
        self.msa_index = msa_index

    def set_amplicon_number(self, amplicon_number) -> None:
        self.amplicon_number = amplicon_number

    def set_pool_number(self, pool_number) -> None:
        self.amplicon_number = pool_number

    def start(self) -> int:
        return min(self.fprimer.starts())

    def end(self) -> int:
        return max(self.rprimer.ends())

    def inter_free(self, cfg) -> bool:
        """
        True means interaction
        """
        return all_inter_checker(self.fprimer.seqs, self.rprimer.seqs, cfg=cfg)
 
    def all_seqs(self) -> set[str]:
        return [x for x in self.fprimer.seqs] + [x for x in self.rprimer.seqs]

    def __hash__(self) -> int:
        return hash(f"{self.start}{self.end}{self.all_seqs()}")

    def __str__(self, ref_name, amplicon_prefix):
        return self.fprimer.__str__(
            referance=f"{ref_name}",
            amplicon_prefix=f"{amplicon_prefix}_{self.amplicon_number}",
            pool=self.pool + 1,
        ) + self.rprimer.__str__(
            referance=f"{ref_name}",
            amplicon_prefix=f"{amplicon_prefix}_{self.amplicon_number}",
            pool=self.pool + 1,
        )
    

class Scheme():
    _pools : list[list[PrimerPair]]
    _current_pool: int
    npools: int
    _last_pp_added: PrimerPair
    cfg: dict

    def __init__(self,cfg):
        self.n_pools = cfg["npools"]
        self._pools = [[] for _ in range(self.n_pools)]
        self._current_pool = 0
        self._pp_number = 1
        self.cfg = cfg
        
    
    @property
    def npools(self) -> int:
        return self.n_pools
    
    def next_pool(self) -> int:
        return (self._current_pool + 1) % self.n_pools
    
    def add_primer_pair_to_pool(self, primerpair, pool):
        primerpair.pool = pool
        primerpair.amplicon_number = len([pp for sublist in self._pools for pp in sublist])
        
        # Adds the primerpair to the spesified pool
        self._pools[pool].append(primerpair)
        self._current_pool = self.next_pool()
        self._last_pp_added = primerpair

    def add_primer_pair(self, primerpair):
        "Adds primerpair to the current pool, and updates the current pool"
        # Adds the primerpair to the current pool
        self.add_primer_pair_to_pool(primerpair, self._current_pool)

    def get_seqs_in_pool(self) -> list[str]:
        return [y for sublist in (x.all_seqs() for x in self._pools[self._current_pool]) for y in sublist]
    
    def get_leading_coverage_edge(self) -> tuple[int,int]:
        """This will return the furthest primer-trimmed region with coverage"""
        # This will crash if no primer has been added, but should not be called until one is
        return self._last_pp_added.rprimer.start
    
    def get_leading_amplicon_edge(self) -> tuple[int,int]:
        """This will return the furthest point of an amplicon"""
        # This will crash if no primer has been added, but should not be called until one is
        return max(self._last_pp_added.rprimer.ends())
    
    def try_ol_primerpairs(self, all_pp_list, thermo_cfg) -> bool:
        """
        This will try and add this primerpair into any valid pool.
        Will return true if the primerpair has been added
        """
        last_pool = self._last_pp_added.pool
        # Find what other pools to look in
        pos_pools_indexes = [(last_pool + i)%self.n_pools for i in range(self.n_pools ) if (last_pool + i)%self.n_pools  != last_pool]

        # Create a hashmap of all sequences in each pool for quick look up
        index_to_seqs: dict[int:list[str]] = {index: [y for sublist in (x.all_seqs() for x in self._pools[index]) for y in sublist] for index in pos_pools_indexes}
        
        # Find pp that could ol, depending on which pool
        pos_ol_pp = [pp for pp in all_pp_list if min(pp.fprimer.starts()) < self.get_leading_coverage_edge() - self.cfg["min_overlap"] and pp.rprimer.start > self.get_leading_amplicon_edge() + self.cfg["min_overlap"]]
        
        #pos_ol_pp = get_pp_window(all_pp_list, fp_start=self.get_leading_coverage_edge() - self.cfg["min_overlap"] - self.cfg["amplicon_size_max"], fp_end=self.get_leading_coverage_edge() - self.cfg["min_overlap"], rp_end=self.get_leading_amplicon_edge() + self.cfg["min_overlap"])
        # Sort the primerpairs depending on how good they are
        ## TODO produce a better score function
        pos_ol_pp.sort(key = lambda pp: (-pp.rprimer.start,len(pp.all_seqs())))

        # For each primerpair
        for ol_pp in pos_ol_pp:
            # For each pool
            for pool_index in pos_pools_indexes:
                # If the pool is empty 
                if not self._pools[pool_index]:
                    self.add_primer_pair_to_pool(ol_pp, pool_index)
                    return True
                # Ensure the new Primerpair doesn't clash with current primers in the pool, either sterically or an interaction  
                if max(self._pools[pool_index][-1].rprimer.ends()) < min(ol_pp.fprimer.starts()) and not all_inter_checker(ol_pp.all_seqs(), index_to_seqs.get(pool_index), thermo_cfg):
                    self.add_primer_pair_to_pool(ol_pp, pool_index)
                    return True

        # If non of the primers work, return false
        return False
    
    def try_walk_primerpair(self, all_pp_list, thermo_cfg) -> bool:
        """
        Find the next valid primerpair while walking forwards
        """
        last_pool = self._last_pp_added.pool
        # Find what other pools to look in, can look in same pool
        pos_pools_indexes = [(last_pool + i)%self.n_pools for i in range(self.n_pools )]

        # Create a hashmap of all sequences in each pool for quick look up
        index_to_seqs: dict[int:list[str]] = {index: [y for sublist in (x.all_seqs() for x in self._pools[index]) for y in sublist] for index in pos_pools_indexes}
        
        # Find all posiable valid primerpairs
        pos_walk_pp = [pp for pp in all_pp_list if pp.fprimer.end > (self.get_leading_coverage_edge() - (self.cfg["min_overlap"] * 2))]
        # Sort walk primers by increasing start position
        # TODO write a better scoring funciton
        pos_walk_pp.sort(key = lambda pp: (pp.fprimer.end,len(pp.all_seqs())))

        # For each primer, try each pool
        for walk_pp in pos_walk_pp:
            for pool_index in pos_pools_indexes:
                # If the pool is empty add the first primer
                if not self._pools[pool_index]:
                    self.add_primer_pair_to_pool(walk_pp, pool_index)
                    return True
                # Check if the walking primer clashes with the primer already in the pool
                if max(self._pools[pool_index][-1].rprimer.ends()) < min(walk_pp.fprimer.starts()) and not all_inter_checker(walk_pp.all_seqs(), index_to_seqs.get(pool_index), thermo_cfg): 
                    self.add_primer_pair_to_pool(walk_pp, pool_index)
                    return True
        
        return False

    def all_primers(self) -> list[PrimerPair]:
        return [pp for pool in (x for x in self._pools) for pp in pool]


    


    