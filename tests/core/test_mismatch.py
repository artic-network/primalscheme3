import unittest

from primalscheme3.core.config import Config
from primalscheme3.core.mismatches import (
    MatchDB,
    detect_new_products,
    detect_products,
    generate_single_mismatches,
)


class Test_DetectNewProducts(unittest.TestCase):
    def test_detect_new_products_inter(self):
        """
        Forms an interaction
        """
        result = detect_new_products({(0, 100, "+")}, {(0, 100, "+"), (0, 200, "-")})
        self.assertTrue(result)

    def test_detect_new_products_no_inter(self):
        """
        No interaction
        """
        result = detect_new_products({(0, 100, "+")}, {(0, 100, "+"), (0, 99, "-")})
        self.assertFalse(result)

    def test_detect_new_products_empty(self):
        """
        Empty
        """
        result = detect_new_products({}, {})  # type: ignore
        self.assertFalse(result)

    def test_detect_new_products_altmsa(self):
        """
        No interaction as the matches are from different MSAs
        """
        result = detect_new_products({(0, 100, "+")}, {(0, 100, "+"), (1, 200, "-")})
        self.assertFalse(result)


class Test_DetectProducts(unittest.TestCase):
    def test_detect_products_inter(self):
        """
        Forms an interaction
        """
        result = detect_products({(0, 100, "+"), (0, 200, "-")})
        self.assertTrue(result)

    def test_detect_products_no_inter(self):
        """
        No interaction
        """
        result = detect_products({(0, 100, "+")})
        self.assertFalse(result)

    def test_detect_products_empty(self):
        """
        Empty
        """
        result = detect_products({})  # type: ignore
        self.assertFalse(result)

    def test_detect_products_altmsa(self):
        """
        No interaction as the matches are from different MSAs
        """
        result = detect_products({(0, 100, "+"), (1, 200, "-")})
        self.assertFalse(result)


def hamming_distance(seq1: str, seq2: str) -> int:
    """
    Calculates the hamming distance between two sequences
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be the same length")
    return sum(base1 != base2 for base1, base2 in zip(seq1, seq2, strict=False))


class Test_MatchDB(unittest.TestCase):
    # Notes for this test
    ## "GTTCAGTATCGACGCGACAA" is at (433, +) in ./tests/test_mismatch.fasta
    ## "TTGTCGCGTCGATACTGAAC" is at (433, -) in ./tests/test_mismatch.fasta

    ## "CTAGCACACTTAAGACGGAG" is found at [(720, +), (780, +)]

    ## 'GCAGGG---TACACTCGGACTCA' -> 'GCAGGGTACACTCGGACTCA' is found in sequence 3 at (496,+)

    ## GGG---TACACTCGGACTCAGGC -> GCCTGAGTCCGAGTGTACCC is found in sequence 3 at (499,-)
    def setUp(self):
        self.matchdb_path = "./tests/core/testcase"
        self.config = Config()

    def test_matchdb_createdb_file(self):
        """
        This tests the creation of the matchdb and the find_match function
        """
        self.config.mismatch_kmersize = 20
        matchdb = MatchDB(
            self.matchdb_path, ["tests/core/test_mismatch.fasta"], self.config
        )

        result = matchdb.find_match("CTAGCACACTTAAGACGGAG")
        self.assertEqual(result, [[0, 720, "+"], [0, 780, "+"]])

    def test_matchdb_find_rc(self):
        """
        This creates a db and then finds the reverse complement of a sequence
        """
        self.config.mismatch_kmersize = 20
        matchdb = MatchDB(
            self.matchdb_path, ["tests/core/test_mismatch.fasta"], self.config
        )

        result = matchdb.find_match("TTGTCGCGTCGATACTGAAC")
        self.assertEqual(result, [[0, 433, "-"]])

    def test_matchdb_find_gap(self):
        """
        This creates a db and then finds a sequence with a gap
        """
        self.config.mismatch_kmersize = 20
        matchdb = MatchDB(
            self.matchdb_path, ["tests/core/test_mismatch.fasta"], self.config
        )

        result = matchdb.find_match("GCAGGGTACACTCGGACTCA")
        self.assertEqual(result, [[0, 496, "+"]])

    def test_matchdb_find_gap_rc(self):
        """
        This creates a db and then finds a sequence with a gap in the reverse complement
        """
        self.config.mismatch_kmersize = 20
        matchdb = MatchDB(
            self.matchdb_path, ["tests/core/test_mismatch.fasta"], self.config
        )

        result = matchdb.find_match("GCCTGAGTCCGAGTGTACCC")
        self.assertEqual(result, [[0, 499, "-"]])


class Test_GenerateSingleMismatches(unittest.TestCase):
    def test_generate_singleton(self):
        """This should return all bases"""
        result = generate_single_mismatches("A")
        self.assertEqual(result, {"A", "C", "G", "T"})

    def test_generate_kmers_invalid_base(self):
        """This should raise an ValueError"""
        self.assertRaises(
            ValueError,
            generate_single_mismatches,
            base_seq="ACTGN",
        )

    def test_generate_kmers_correct_distance(self):
        """All seqs should have distance of 1 or 0"""
        input_seq = "ATGCATGATCGAC"
        single_mismatches = generate_single_mismatches(input_seq)

        # Generate the hamming distances for all seqs
        hamming_distances = {
            hamming_distance(input_seq, seq) for seq in single_mismatches
        }

        self.assertEqual(hamming_distances, {0, 1})


if __name__ == "__main__":
    unittest.main()
