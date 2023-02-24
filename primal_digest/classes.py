from primal_digest.config import AMBIGUOUS_DNA_COMPLEMENT
from primal_digest.iteraction import all_inter_checker

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