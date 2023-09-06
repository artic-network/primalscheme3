import unittest
from primal_digest.digestion import (
    reduce_kmers,
    hamming_dist,
    walk_left,
    walk_right,
    wrap_walk,
    mp_r_digest,
    mp_f_digest,
    digest,
)
from primal_digest.classes import FKmer, RKmer, PrimerPair
from primal_digest.config import config_dict as cfg
from primal_digest.thermo import calc_tm

from primal_digest.errors import *
import numpy as np


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


class Test_WalkLeft(unittest.TestCase):
    def test_walk_left(self):
        """
        Ensure the Tm prodvided is greater than the min
        """
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]
        array_list = []
        for seq in seqs:
            array_list.append([x for x in seq])

        msa_array = np.array(array_list)
        right = 60
        left = 60 - 15

        result = wrap_walk(
            walk_left,
            msa_array,
            right,
            left,
            0,
            "".join(msa_array[0, left:right]),
            cfg,
        )

        expected = "GTCCAATGGTGCAAAAGGTATAATCATTAAT"

        self.assertGreaterEqual(
            calc_tm(expected, cfg),
            calc_tm([x for x in result][0], cfg),
        )

    def test_walk_left_expandambs(self):
        """
        Tests if ambs are expanded correctly
        """
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTAYAATCATTAATGT",
        ]
        array_list = []
        for seq in seqs:
            array_list.append([x for x in seq])

        msa_array = np.array(array_list)
        right = 60
        left = 60 - 15
        result = walk_left(
            msa_array,
            right,
            left,
            0,
            "".join(msa_array[0, left:right]),
            cfg,
        )

        expected = ["CCAATGGTGCAAAAGGTACAATCATTAAT", "GTCCAATGGTGCAAAAGGTATAATCATTAAT"]

        # Ensure both a sorted
        result = [x for x in result]
        result.sort()
        expected.sort()
        self.assertEqual(result, expected)

    def test_walk_left_containsinvalid(self):
        """
        Expect to return None of non valid bases
        """
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTA..AATCATTAATGT",
        ]
        array_list = []
        for seq in seqs:
            array_list.append([x for x in seq])

        msa_array = np.array(array_list)
        right = 60
        left = 60 - 15

        result = wrap_walk(
            walk_left,
            msa_array,
            right,
            left,
            0,
            "".join(msa_array[0, left:right]),
            cfg,
        )

        expected = [ContainsInvalidBase()]
        self.assertEqual(result, expected)

    def test_walk_left_outsidearray(self):
        """
        Test None is returned if the index walks outside the array
        """
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]
        array_list = []
        for seq in seqs:
            array_list.append([x for x in seq])

        msa_array = np.array(array_list)
        right = 15
        left = right - 10
        result = wrap_walk(
            walk_left,
            msa_array,
            right,
            left,
            0,
            "".join(msa_array[0, left:right]),
            cfg,
        )

        expected = [WalksOut()]
        self.assertEqual(result, expected)


class Test_WalkRight(unittest.TestCase):
    def test_walk_right(self):
        """
        Ensure the Tm prodvided is greater than the min
        """
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]
        array_list = []
        for seq in seqs:
            array_list.append([x for x in seq])

        msa_array = np.array(array_list)
        left = 10
        right = left + 10

        result = wrap_walk(
            walk_right,
            msa_array,
            right,
            left,
            0,
            "".join(msa_array[0, left:right]),
            cfg,
        )

        expected = "GTCCAATGGTGCAAAAGGTATAATCATTAAT"

        self.assertGreaterEqual(
            calc_tm(expected, cfg),
            calc_tm([x for x in result][0], cfg),
        )

    def test_walk_right_expandambs(self):
        """
        Tests if ambs are expanded correctly
        """
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTAYAATCATTAATGT",
        ]
        array_list = []
        for seq in seqs:
            array_list.append([x for x in seq])

        msa_array = np.array(array_list)
        left = 30
        right = left + 10

        result = wrap_walk(
            walk_right,
            msa_array,
            right,
            left,
            0,
            "".join(msa_array[0, left:right]),
            cfg,
        )

        expected = ["TCCAATGGTGCAAAAGGTACAATCA", "TCCAATGGTGCAAAAGGTATAATCATTAATG"]

        # Ensure both a sorted
        result = [x for x in result]
        result.sort()
        expected.sort()
        self.assertEqual(result, expected)

    def test_walk_left_containsinvalid(self):
        """
        Expect to return None of non valid bases
        """
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTA..AATCATTAATGT",
        ]
        array_list = []
        for seq in seqs:
            array_list.append([x for x in seq])

        msa_array = np.array(array_list)
        left = 30
        right = left + 10
        result = wrap_walk(
            walk_right,
            msa_array,
            right,
            left,
            0,
            "".join(msa_array[0, left:right]),
            cfg,
        )

        expected = [ContainsInvalidBase()]
        self.assertEqual(result, expected)

    def test_walk_left_contains1invalid(self):
        """
        Expect to return None of non valid bases
        """
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTA..AATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTAAAAATCATTAATGT",
        ]
        array_list = []
        for seq in seqs:
            array_list.append([x for x in seq])

        msa_array = np.array(array_list)
        left = 30
        right = left + 10

        result = wrap_walk(
            walk_right,
            msa_array,
            right,
            left,
            0,
            "".join(msa_array[0, left:right]),
            cfg,
        )

        # Process the second sequence in a hacky way
        result.extend(
            wrap_walk(
                walk_right,
                msa_array,
                right,
                left,
                1,
                "".join(msa_array[0, left:right]),
                cfg,
            )
        )

        expected = [ContainsInvalidBase(), "TCCAATGGTGCAAAAGGTAAAAATCATTAA"]
        self.assertEqual(result, expected)

    def test_walk_left_outsidearray(self):
        """
        Test None is returned if the index walks outside the array
        """
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]
        array_list = []
        for seq in seqs:
            array_list.append([x for x in seq])

        msa_array = np.array(array_list)
        left = 50
        right = left + 10
        result = wrap_walk(
            walk_right,
            msa_array,
            right,
            left,
            0,
            "".join(msa_array[0, left:right]),
            cfg,
        )

        expected = [WalksOut()]
        self.assertEqual(result, expected)


class Test_MPRDigest(unittest.TestCase):
    def test_mp_r_digest(self):
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]
        array_list = []
        for seq in seqs:
            array_list.append([x for x in seq])

        msa_array = np.array(array_list)
        cfg["reducekmers"] = False
        cfg["dimerscore"] = -26

        data = (msa_array, cfg, 20, 0)
        result = mp_r_digest(data)

        # The Expected Sequence
        expected = {"ACCTTTTGCACCATTGGACATTAATGAT"}

        self.assertEqual(result.seqs, expected)

    def test_mp_r_digest_one_invalid(self):
        """The invalid base and min_freq 0 should return None"""
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCANTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]
        array_list = []
        for seq in seqs:
            array_list.append([x for x in seq])

        msa_array = np.array(array_list)
        cfg["reducekmers"] = False
        cfg["dimerscore"] = -26

        data = (msa_array, cfg, 20, 0)
        result = mp_r_digest(data)

        # The Expected Sequence
        expected = None

        self.assertEqual(result, expected)

    def test_mp_r_digest_one_invalid(self):
        """The invalid base and min_freq 0.5 should return sequence"""
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCANTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]
        array_list = []
        for seq in seqs:
            array_list.append([x for x in seq])

        msa_array = np.array(array_list)
        cfg["reducekmers"] = False
        cfg["dimerscore"] = -26

        data = (msa_array, cfg, 20, 0.5)
        result = mp_r_digest(data)

        # The Expected Sequence
        expected = {"ACCTTTTGCACCATTGGACATTAATGAT"}

        self.assertEqual(result.seqs, expected)


class Test_MPFDigest(unittest.TestCase):
    def test_mp_f_digest(self):
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]
        array_list = []
        for seq in seqs:
            array_list.append([x for x in seq])

        msa_array = np.array(array_list)
        cfg["reducekmers"] = False
        cfg["dimerscore"] = -26

        data = (msa_array, cfg, 40, 0)
        result = mp_f_digest(data)

        # The Expected Sequence
        expected = {"GCAAAAGGTATAATCATTAATGTCCAATGGTG"}

        self.assertEqual(result.seqs, expected)

    def test_mp_f_digest_one_invalid(self):
        """The invalid base and min_freq 0 should return None"""
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCANTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]
        array_list = []
        for seq in seqs:
            array_list.append([x for x in seq])

        msa_array = np.array(array_list)
        cfg["reducekmers"] = False

        data = (msa_array, cfg, 40, 0)
        result = mp_f_digest(data)

        # The Expected Sequence
        expected = None

        self.assertEqual(result, expected)

    def test_mp_r_digest_one_invalid(self):
        """The invalid base and min_freq 0.5 should return sequence"""
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCANTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]
        array_list = []
        for seq in seqs:
            array_list.append([x for x in seq])

        msa_array = np.array(array_list)
        cfg["reducekmers"] = False
        cfg["dimerscore"] = -26

        data = (msa_array, cfg, 40, 0.5)
        result = mp_f_digest(data)

        # The Expected Sequence
        expected = {"GCAAAAGGTATAATCATTAATGTCCAATGGTG"}

        self.assertEqual(result.seqs, expected)


class TestDigest(unittest.TestCase):
    def test_digestion_valid_fkmer(self):
        """Test the digestion"""
        self.maxDiff = None
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCANTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]
        array_list = []
        for seq in seqs:
            array_list.append([x for x in seq])

        msa_array = np.array(array_list)
        cfg["reducekmers"] = False
        cfg["minbasefreq"] = 0
        cfg["dimerscore"] = -26

        results = digest(msa_array=msa_array, cfg=cfg, indexes=([60], []))
        expected_fkmer = FKmer(60, {"GTCCAATGGTGCAAAAGGTATAATCATTAAT"})

        fkmer = results[0][0]
        self.assertEqual(fkmer, expected_fkmer)

    def test_digestion_valid_rkmer(self):
        """Test the digestion"""
        self.maxDiff = None
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCANTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]
        array_list = []
        for seq in seqs:
            array_list.append([x for x in seq])

        msa_array = np.array(array_list)
        cfg["reducekmers"] = False
        cfg["minbasefreq"] = 0
        cfg["dimerscore"] = -26

        results = digest(msa_array=msa_array, cfg=cfg, indexes=([], [25]))
        expected_rkmer = RKmer(25, {'TGATTATACCTTTTGCACCATTGGACATTA'})

        rkmer = results[1][0]
        print(rkmer.start)
        print(rkmer.seqs)
        self.assertEqual(rkmer, expected_rkmer)


if __name__ == "__main__":
    unittest.main()
