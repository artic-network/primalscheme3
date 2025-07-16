# Module imports
from primalschemers._core import FKmer, RKmer, do_pool_interact  # type: ignore

from primalscheme3.core.config import Config
from primalscheme3.core.thermo import calc_tm, gc


class PrimerPair:
    fprimer: FKmer
    rprimer: RKmer
    amplicon_number: int
    pool: int
    msa_index: int
    chrom_name: str | None
    amplicon_prefix: str | None
    _score: float | None

    __slots__ = [
        "fprimer",
        "rprimer",
        "amplicon_number",
        "pool",
        "msa_index",
        "chrom_name",
        "amplicon_prefix",
        "_score",
    ]

    def __init__(
        self,
        fprimer,
        rprimer,
        msa_index,
        amplicon_number=-1,
        pool=-1,
    ):
        self.fprimer = fprimer
        self.rprimer = rprimer
        self.amplicon_number = amplicon_number
        self.pool = pool
        self.msa_index = msa_index
        self.chrom_name = None
        self.amplicon_prefix = None
        self._score = None

    def get_score(self, target_gc=0.5):
        """
        Returns the mean gc diff of the primerpair
        """
        if self._score is None:
            self._score = sum(
                [abs(target_gc - (gc(x) / 100)) for x in self.all_seqs()]
            ) / len(self.all_seqs())
        return self._score

    def regions(self) -> tuple[tuple[int, int], tuple[int, int]]:
        return self.fprimer.region(), self.rprimer.region()

    def set_amplicon_number(self, amplicon_number) -> None:
        self.amplicon_number = amplicon_number

    def set_pool_number(self, pool_number) -> None:
        self.amplicon_number = pool_number

    def find_matches(self, matchDB, fuzzy, remove_expected, kmersize) -> set[tuple]:
        """
        Find matches for the FKmer and RKmer
        """
        matches = set()
        # Find the FKmer matches
        matches.update(
            matchDB.find_fkmer(
                self.fprimer,
                fuzzy=fuzzy,
                remove_expected=remove_expected,
                kmersize=kmersize,
                msaindex=self.msa_index,
            )
        )
        # Find the RKmer matches
        matches.update(
            matchDB.find_rkmer(
                self.rprimer,
                fuzzy=fuzzy,
                remove_expected=remove_expected,
                kmersize=kmersize,
                msaindex=self.msa_index,
            )
        )
        return matches

    def kmers(self):
        """
        Returns the FKmer and RKmer
        """
        return self.fprimer, self.rprimer

    def primertrimmed_region(self) -> tuple[int, int]:
        """
        Returns the region of the primertrimed region
        Right position is non-inclusive
        """
        return self.fprimer.end, self.rprimer.start

    def inter_free(self, cfg) -> bool:
        """
        True means interaction
        """
        return do_pool_interact(
            [*self.fprimer.seqs], [*self.rprimer.seqs], cfg["dimerscore"]
        )

    def all_seqs(self) -> list[str]:
        return self.fprimer.seqs() + self.rprimer.seqs()

    def all_seq_bytes(self) -> list[bytes]:
        return self.fprimer.seqs_bytes() + self.rprimer.seqs_bytes()

    def calc_tm(self, config: Config) -> list[float]:
        """
        Calculates the tm for all primers in the PrimerPair
        :param cfg: config dict
        :return: list of tm values
        """
        return [
            calc_tm(
                seq,
                mv_conc=config.mv_conc,
                dv_conc=config.dv_conc,
                dna_conc=config.dna_conc,
                dntp_conc=config.dna_conc,
            )
            for seq in self.all_seqs()
        ]

    def __hash__(self) -> int:
        return hash(f"{self.regions()}{self.all_seqs()}")

    def __eq__(self, other):
        if isinstance(other, PrimerPair):
            return self.__hash__() == other.__hash__()
        else:
            return False

    def to_bed(self) -> str:
        """
        Turns the primerpair into a string for a bed file
        :param chromname: name of the chromosome
        :param amplicon_prefix: prefix for the amplicon
        :return: string for the bed file
        """
        return self.__str__()

    def __str__(self):
        return self.fprimer.to_bed(
            chrom=f"{self.chrom_name}",
            amplicon_prefix=f"{self.amplicon_prefix}_{self.amplicon_number}",
            pool=self.pool + 1,
        ) + self.rprimer.to_bed(
            chrom=f"{self.chrom_name}",
            amplicon_prefix=f"{self.amplicon_prefix}_{self.amplicon_number}",
            pool=self.pool + 1,
        )
