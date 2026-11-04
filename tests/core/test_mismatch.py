import gzip
import pathlib
import random
import tempfile
import unittest

from primalschemers import FKmer, PoolHits, RKmer  # type: ignore

from primalscheme3.core.classes import PrimerPair
from primalscheme3.core.config import Config, MappingType
from primalscheme3.core.mismatches import MatchDB, expected_site
from primalscheme3.core.seq_functions import reverse_complement

FIXTURE = "tests/core/test_mismatch.fasta"


def random_seq(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


class TestPoolHits(unittest.TestCase):
    def test_new_plus_with_pool_minus(self):
        pool = PoolHits()
        pool.add([(0, 0, 100, "+"), (0, 0, 200, "-")])
        self.assertTrue(pool.forms_product([(0, 0, 100, "+")], 2000))

    def test_minus_upstream_of_plus(self):
        pool = PoolHits()
        pool.add([(0, 0, 100, "+"), (0, 0, 99, "-")])
        self.assertFalse(pool.forms_product([(0, 0, 100, "+")], 2000))

    def test_empty(self):
        self.assertFalse(PoolHits().forms_product([], 2000))

    def test_other_reference(self):
        pool = PoolHits()
        pool.add([(1, 0, 200, "-")])
        self.assertFalse(pool.forms_product([(0, 0, 100, "+")], 2000))

    def test_product_size(self):
        pool = PoolHits()
        pool.add([(0, 0, 200, "-")])
        self.assertTrue(pool.forms_product([(0, 0, 100, "+")], 101))
        self.assertFalse(pool.forms_product([(0, 0, 100, "+")], 100))

    def test_new_hits_with_each_other(self):
        f, r = (0, 0, 100, "+"), (0, 0, 200, "-")
        self.assertTrue(PoolHits().forms_product([f, r], 2000))
        self.assertFalse(PoolHits().forms_product([f, r], 2000, [f, r]))

    def test_remove_is_multiset(self):
        pool = PoolHits()
        pool.add([(0, 0, 100, "+")])
        pool.add([(0, 0, 100, "+")])
        pool.remove([(0, 0, 100, "+")])
        self.assertEqual(len(pool), 1)
        self.assertTrue(pool.forms_product([(0, 0, 200, "-")], 2000))

    def test_invalid_strand(self):
        with self.assertRaises(ValueError):
            PoolHits().add([(0, 0, 100, "x")])


class TestMatchDBFixture(unittest.TestCase):
    """Positions are columns (consensus mapping) of the primer's 3' end base"""

    def setUp(self):
        self.config = Config(mapping=MappingType.CONSENSUS)
        self.config.mismatch_kmersize = 20

    def test_describe_unused(self):
        self.assertEqual(
            MatchDB([], [], self.config).describe(), "Mismatch database not used"
        )

    def test_find_plus(self):
        matchdb = MatchDB([FIXTURE], [], self.config)
        kmer = FKmer([b"CTAGCACACTTAAGACGGAG"], 740)
        self.assertEqual(matchdb.find_kmer(kmer), [(0, 0, 739, "+"), (0, 0, 799, "+")])

    def test_find_minus(self):
        matchdb = MatchDB([FIXTURE], [], self.config)
        kmer = RKmer([b"TTGTCGCGTCGATACTGAAC"], 433)
        self.assertEqual(matchdb.find_kmer(kmer), [(0, 0, 433, "-")])
        self.assertEqual(expected_site(kmer, 0), (0, 0, 433, "-"))

    def test_find_across_gap(self):
        # Sequence 3 has GCAGGG---TACACTCGGACTCA from column 496
        matchdb = MatchDB([FIXTURE], [], self.config)
        self.assertEqual(
            matchdb.find_kmer(FKmer([b"GCAGGGTACACTCGGACTCA"], 0)),
            [(0, 0, 518, "+")],
        )
        self.assertEqual(
            matchdb.find_kmer(RKmer([b"GCCTGAGTCCGAGTGTACCC"], 0)),
            [(0, 0, 499, "-")],
        )

    def test_fuzzy(self):
        mutated = b"CTAGCACACTTAAGACGGAC"
        self.config.mismatch_fuzzy = False
        self.assertEqual(
            MatchDB([FIXTURE], [], self.config).find_kmer(FKmer([mutated], 0)), []
        )
        self.config.mismatch_fuzzy = True
        self.assertEqual(
            MatchDB([FIXTURE], [], self.config).find_kmer(FKmer([mutated], 0)),
            [(0, 0, 739, "+"), (0, 0, 799, "+")],
        )


class TestMatchDBSynthetic(unittest.TestCase):
    """A random genome as the MSA, and background genomes built from it"""

    def setUp(self):
        rng = random.Random(42)
        self.genome = random_seq(rng, 1000)
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = pathlib.Path(self.tmp.name)
        self.msa = self.dir / "msa.fasta"
        self.msa.write_text(f">genome\n{self.genome}\n")

        # A pair whose primers are longer than the kmer: only 3' ends count
        self.fkmer = FKmer([self.genome[100:125].encode()], 125)
        self.rkmer = RKmer([reverse_complement(self.genome[400:425]).encode()], 400)
        self.pair = PrimerPair(self.fkmer, self.rkmer, 0)

        # A background carrying both primer sites 200 bases apart, on contig 2
        self.background = self.dir / "background.fasta.gz"
        decoy = (
            random_seq(rng, 50)
            + self.genome[100:125]
            + random_seq(rng, 200)
            + self.genome[400:425]
            + random_seq(rng, 50)
        )
        with gzip.open(self.background, "wt") as f:
            f.write(f">unrelated\n{random_seq(rng, 500)}\n>decoy x\n{decoy}\n")

        self.config = Config()

    def tearDown(self):
        self.tmp.cleanup()

    def test_expected_sites(self):
        matchdb = MatchDB([self.msa], [], self.config)
        self.assertEqual(matchdb.find_kmer(self.fkmer), [(0, 0, 124, "+")])
        self.assertEqual(matchdb.find_kmer(self.rkmer), [(0, 0, 400, "-")])
        self.assertEqual(
            self.pair.expected_matches(), [(0, 0, 124, "+"), (0, 0, 400, "-")]
        )

    def test_reverse_primer_uses_3_prime_end(self):
        # Changing the reverse primer's 5' end doesn't change its hits
        seq = reverse_complement(self.genome[400:425])
        rkmer = RKmer([("AAAAAA" + seq[6:]).encode()], 400)
        matchdb = MatchDB([self.msa], [], self.config)
        self.assertEqual(matchdb.find_kmer(rkmer), [(0, 0, 400, "-")])

    def test_pair_alone_is_not_a_product(self):
        matchdb = MatchDB([self.msa], [], self.config)
        self.assertFalse(
            PoolHits().forms_product(
                self.pair.find_matches(matchdb),
                self.config.mismatch_product_size,
                self.pair.expected_matches(),
            )
        )

    def test_background_product(self):
        matchdb = MatchDB([self.msa], [str(self.background)], self.config)
        self.assertEqual(matchdb.prescreen_msas({}, [self.pair]), 2)
        hits = self.pair.find_matches(matchdb)
        # Decoy is contig 1 of ref 1: F 3' end at 50 + 24, R at 50 + 25 + 200
        self.assertIn((1, 1, 74, "+"), hits)
        self.assertIn((1, 1, 275, "-"), hits)
        self.assertTrue(
            PoolHits().forms_product(
                hits,
                self.config.mismatch_product_size,
                self.pair.expected_matches(),
            )
        )

    def test_cross_msa_hits_kept(self):
        # The same genome as a second MSA: the pair's sites there are off-target
        matchdb = MatchDB([self.msa, self.msa], [], self.config)
        hits = self.pair.find_matches(matchdb)
        self.assertIn((1, 0, 124, "+"), hits)
        self.assertTrue(
            PoolHits().forms_product(
                hits,
                self.config.mismatch_product_size,
                self.pair.expected_matches(),
            )
        )

    def test_missing_background(self):
        with self.assertRaises(ValueError):
            MatchDB([self.msa], [str(self.dir / "missing.fasta")], self.config)


class TestConfig(unittest.TestCase):
    def test_product_size_defaults_to_amplicon_size_max(self):
        config = Config(amplicon_size=400)
        self.assertEqual(config.mismatch_product_size, config.amplicon_size_max)


if __name__ == "__main__":
    unittest.main()
