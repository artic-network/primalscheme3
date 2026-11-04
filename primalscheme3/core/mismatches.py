import pathlib
from collections.abc import Iterable

from primalschemers import FKmer, MisprimeDB, RKmer  # type: ignore

from primalscheme3.core.config import Config, MappingType

# (ref_id, contig_id, pos, strand): pos is the primer's 3' end base
Hit = tuple[int, int, int, str]


class MatchDB:
    """
    Where primers' 3' terminal kmers (config.mismatch_kmersize bases, with one
    mismatch if config.mismatch_fuzzy) bind in the MSAs and background genomes.

    MSAs are refs 0..n_msas, one contig each, with positions in the kmers'
    coordinates. Background genomes follow, one contig per record.
    """

    def __init__(
        self,
        msa_paths: list[str],
        background_paths: list[str],
        config: Config,
    ) -> None:
        self.n_msas = len(msa_paths)
        self.n_backgrounds = len(background_paths)
        self._db = MisprimeDB(
            [str(p) for p in msa_paths],
            [str(p) for p in background_paths],
            config.mismatch_kmersize,
            config.mismatch_fuzzy,
            config.mapping == MappingType.FIRST,
        )

    def describe(self) -> str:
        """What is screened, for the log"""
        if self.n_msas == 0 and self.n_backgrounds == 0:
            return "Mismatch database not used"
        return (
            f"Mismatch database: {self.n_msas} MSA(s), "
            f"{self.n_backgrounds} background genome(s)"
        )

    def prescreen(self, kmers: Iterable[FKmer | RKmer]) -> int:
        """
        Screen every kmer in one pass over the references, so later lookups
        are cached. Returns how many terminal kmers were screened.
        """
        return self._db.screen([s for kmer in kmers for s in kmer.seqs_bytes()])

    def prescreen_msas(self, msa_dict: dict, extra_pairs: Iterable = ()) -> int:
        """
        Screen every digested primer in the MSAs, and the extra
        primer pairs (e.g. from an input bedfile)
        """
        kmers = [k for msa in msa_dict.values() for k in (*msa.fkmers, *msa.rkmers)]
        kmers += [k for pp in extra_pairs for k in pp.mispriming_kmers()]
        return self.prescreen(kmers)

    def ref_names(self) -> list[str]:
        """Each reference's file name, MSAs first"""
        return [pathlib.Path(p).name for p in self._db.ref_names()]

    def site_counts(self) -> list[tuple[str, int]]:
        """Each reference's path and its distinct binding sites so far"""
        return list(zip(self._db.ref_names(), self._db.site_counts(), strict=True))

    def find_kmer(self, kmer: FKmer | RKmer) -> list[Hit]:
        """The distinct sites any of the kmer's sequences bind"""
        return sorted({h for hits in self._db.hits(kmer.seqs_bytes()) for h in hits})


def expected_site(kmer: FKmer | RKmer, msa_index: int) -> Hit:
    """The kmer's intended binding site, as a hit"""
    if isinstance(kmer, FKmer):
        return (msa_index, 0, kmer.end - 1, "+")
    return (msa_index, 0, kmer.start, "-")
