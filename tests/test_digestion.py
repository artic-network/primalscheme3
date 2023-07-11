import unittest
from primal_digest.digestion import reduce_kmers, hamming_dist
from primal_digest.classes import FKmer, RKmer, PrimerPair


class Test_HammingDist(unittest.TestCase):
    def test_hamming_dist_exact(self):
        s1 = "ATCATGCTAGCTAGCGTA"
        s2 = "ATCATGCTAGCTAGCGTA"

        expected = 0
        result = hamming_dist(s1, s2)

        self.assertEqual(expected, result)

    def test_hamming_dist_mismatch(self):
        s1 = "ACGATCGATCGTAGCTTAGGCTAC"
        s2 = "ACGATCGACCGTAGCTTAGGCTAC"

        expected = 1
        result = hamming_dist(s1, s2)

        self.assertEqual(expected, result)

    def test_hamming_dist_indel(self):
        # Hamming_dist starts at the left side of the string

        #  ACGATCGATCGTAGCTTAGGCTA
        #  ...............|..|....
        # ACGATCGACCGTAGCTTAGGCTAC
        s1 = "ACGATCGATCGTAGCTTAGGCTA"
        s2 = "ACGATCGACCGTAGCTTAGGCTAC"

        expected = 21
        result = hamming_dist(s1, s2)

        self.assertEqual(expected, result)


class Test_ReduceKmers(unittest.TestCase):
    def test_reduce_kmers(self):
        seqs = {
            "CAATGGTGCGAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGAATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAATGGTATAATCATTAATGT",
            "CCAGTGGTGCAAAAGGTATAATCATTAATGT",
        }

        expected = ["CCAATGGTGCAAAAGGTATAATCATTAATGT"]
        result = reduce_kmers(seqs, 1, 6)

        self.assertEqual(expected, list(result))

    def test_reduce_kmers_single(self):
        seqs = {
            "CCAATGGTGCAAAAGGTATAATCATTAATGT",
        }

        expected = ["CCAATGGTGCAAAAGGTATAATCATTAATGT"]
        result = reduce_kmers(seqs, 1, 6)

        self.assertEqual(expected, list(result))

    def test_reduce_kmers_all_dif(self):
        seqs = {
            "CCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATG",
            "CCAATGGTGCAAAAGGTATAATCATTAAT",
            "CCAATGGTGCAAAAGGTATAATCATTAA",
        }

        expected = list(seqs)
        expected.sort()

        result = list(reduce_kmers(seqs, 1, 6))
        result.sort()

        self.assertEqual(expected, list(result))


if __name__ == "__main__":
    unittest.main()
