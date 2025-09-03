import unittest

import numpy as np
from primalschemers import FKmer, RKmer  # type: ignore

from primalscheme3.core.config import Config
from primalscheme3.core.digestion import (
    DIGESTION_ERROR,
    digest,
    f_digest_index,
    hamming_dist,
    parse_error,
    r_digest_index,
    walk_left,
    walk_right,
    wrap_walk,
)
from primalscheme3.core.errors import (
    ContainsInvalidBase,
    CustomRecursionError,
    GapOnSetBase,
    WalksOut,
    WalksTooFar,
)
from primalscheme3.core.progress_tracker import ProgressManager as PM
from primalscheme3.core.thermo import calc_tm


class Test_parse_errors(unittest.TestCase):
    def test_parse_error(self):
        """
        Test the parsing of errors
        """
        # Test that ContainsInvalidBase() is detected
        self.assertEqual(
            parse_error({"ATCG", ContainsInvalidBase()}),
            DIGESTION_ERROR.CONTAINS_INVALID_BASE,
        )

        # Test that WalksOut() is detected
        self.assertEqual(
            parse_error({"ATCG", WalksOut()}),
            DIGESTION_ERROR.WALKS_OUT,
        )

        # Test that GapOnSetBase() is detected
        self.assertEqual(
            parse_error({"ATCG", GapOnSetBase()}), DIGESTION_ERROR.GAP_ON_SET_BASE
        )

        # CustomRecursionError() is detected
        self.assertEqual(
            parse_error({"ATCG", CustomRecursionError()}),
            DIGESTION_ERROR.CUSTOM_RECURSION_ERROR,
        )

        # WalksOut() is detected
        self.assertEqual(
            parse_error({"ATCG", WalksTooFar()}),
            DIGESTION_ERROR.WALK_TO_FAR,
        )

        # Test that DIGESTION_ERROR.AMB_FAIL is returned on no error provided
        self.assertEqual(parse_error({"ATCG"}), DIGESTION_ERROR.AMB_FAIL)


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


class Test_WalkLeft(unittest.TestCase):
    config = Config()

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
            self.config,
        )

        expected = "GTCCAATGGTGCAAAAGGTATAATCATTAAT"

        self.assertGreaterEqual(
            calc_tm(
                expected,
                mv_conc=self.config.mv_conc,
                dv_conc=self.config.dv_conc,
                dna_conc=self.config.dna_conc,
                dntp_conc=self.config.dntp_conc,
            ),
            calc_tm(
                [x for x in result][0],
                mv_conc=self.config.mv_conc,
                dv_conc=self.config.dv_conc,
                dna_conc=self.config.dna_conc,
                dntp_conc=self.config.dntp_conc,
            ),
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
        right = 65
        left = right - 10
        result = walk_left(
            msa_array,
            right,
            left,
            0,
            "".join(msa_array[0, left:right]),
            self.config,
        )

        expected = ["CAATGGTGCAAAAGGTATAATCATTAATGT", "TGGTGCAAAAGGTACAATCATTAATGT"]

        # Ensure both a sorted
        result = [x for x in result]  # type: ignore
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
        right = 65
        left = right - 10

        result = wrap_walk(
            walk_left,
            msa_array,
            right,
            left,
            0,
            "".join(msa_array[0, left:right]),
            self.config,
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
            self.config,
        )

        expected = [WalksOut()]
        self.assertEqual(result, expected)


class Test_WalkRight(unittest.TestCase):
    config = Config()

    def test_walk_right(self):
        """
        Ensure the Tm provided is greater than the min
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
            self.config,
        )

        expected = "GTCCAATGGTGCAAAAGGTATAATCATTAAT"

        self.assertGreaterEqual(
            calc_tm(
                expected,
                mv_conc=self.config.mv_conc,
                dv_conc=self.config.dv_conc,
                dna_conc=self.config.dna_conc,
                dntp_conc=self.config.dntp_conc,
            ),
            calc_tm(
                [x for x in result][0],
                mv_conc=self.config.mv_conc,
                dv_conc=self.config.dv_conc,
                dna_conc=self.config.dna_conc,
                dntp_conc=self.config.dntp_conc,
            ),
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
            self.config,
        )

        expected = ["TCCAATGGTGCAAAAGGTACAATC", "TCCAATGGTGCAAAAGGTATAATCAT"]

        # Ensure both a sorted
        result = [x for x in result]
        result.sort()  # type: ignore
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
            self.config,
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
            self.config,
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
                self.config,
            )
        )

        expected = [ContainsInvalidBase(), "TCCAATGGTGCAAAAGGTAAAAATCA"]
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
            self.config,
        )

        expected = [WalksOut()]
        self.assertEqual(result, expected)


class Test_MPRDigest(unittest.TestCase):
    def setUp(self):
        self.config: Config = Config()

    def create_array(self, seqs) -> np.ndarray:
        array_list = []
        for seq in seqs:
            array_list.append([x for x in seq])
        return np.array(array_list)

    def test_r_digest_index(self):
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]

        result = r_digest_index(self.create_array(seqs), self.config, 20, 0)

        # The Expected Sequence
        expected = ["CTTTTGCACCATTGGACATTAATGAT"]

        self.assertEqual(result.seqs(), expected)  # type: ignore

    def test_r_digest_index_one_invalid_0(self):
        """The invalid base and min_freq 0 should return None"""
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCANTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]
        result = r_digest_index(self.create_array(seqs), self.config, 20, 0)

        # The Expected Sequence
        expected = (20, DIGESTION_ERROR.CONTAINS_INVALID_BASE)

        self.assertEqual(result, expected)

    def test_r_digest_index_one_invalid_MBF(self):
        """The invalid base and min_freq 0.5 should return sequence"""
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCANTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]

        result = r_digest_index(self.create_array(seqs), self.config, 20, 0.5)

        # The Expected Sequence
        expected = ["CTTTTGCACCATTGGACATTAATGAT"]

        self.assertEqual(result.seqs(), expected)  # type: ignore

    def test_r_digest_index_walkout(self):
        """
        Tests that walking out of the array returns the correct error
        """
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCANTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]
        result = r_digest_index(self.create_array(seqs), self.config, 60, 0.5)

        # The Expected Sequence
        expected = (60, DIGESTION_ERROR.WALKS_OUT)
        self.assertEqual(result, expected)  # type: ignore

    def test_r_digest_index_gaponsetbase(self):
        """
        Tests GapOnSetBase is returned
        """
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCA-TAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCA-TAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]
        result = r_digest_index(self.create_array(seqs), self.config, 24, 0)

        # The Expected Sequence
        expected = (24, DIGESTION_ERROR.GAP_ON_SET_BASE)
        self.assertEqual(result, expected)  # type: ignore

    def test_r_digest_index_walktofar(self):
        """
        Tests WalksToFar is returned
        """

        seqs = [
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]

        local_cfg = Config()
        local_cfg.primer_max_walk = 10  # Force the primer to walk to far

        result = r_digest_index(self.create_array(seqs), local_cfg, 10, 0)

        # The Expected Sequence
        expected = (10, DIGESTION_ERROR.WALK_TO_FAR)
        self.assertEqual(result, expected)  # type: ignore


class Test_MPFDigest(unittest.TestCase):
    def setUp(self):
        self.config = Config()

    def create_array(self, seqs) -> np.ndarray:
        array_list = []
        for seq in seqs:
            array_list.append([x for x in seq])
        return np.array(array_list)

    def test_f_digest_index(self):
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]

        result = f_digest_index(self.create_array(seqs), self.config, 40, 0)

        # The Expected Sequence
        expected = ["AAAAGGTATAATCATTAATGTCCAATGGTG"]

        self.assertEqual(result.seqs(), expected)  # type: ignore

    def test_f_digest_index_one_invalid(self):
        """The invalid base and min_freq 0 should return None"""
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCANTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]
        result = f_digest_index(self.create_array(seqs), self.config, 40, 0)

        # The Expected Sequence
        expected = (40, DIGESTION_ERROR.CONTAINS_INVALID_BASE)

        self.assertEqual(result, expected)

    def test_f_digest_index_one_invalid_with_mbf(self):
        """The invalid base and min_freq 0.5 should return sequence"""
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCANTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]
        result = f_digest_index(self.create_array(seqs), self.config, 40, 0.5)

        # The Expected Sequence
        expected = ["AAAAGGTATAATCATTAATGTCCAATGGTG"]

        self.assertEqual(result.seqs(), expected)  # type: ignore

    def test_f_digest_index_walkout(self):
        """
        Tests that walking out of the array returns the correct error
        """
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCANTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]
        result = f_digest_index(self.create_array(seqs), self.config, 5, 0.5)

        # The Expected Sequence
        expected = (5, DIGESTION_ERROR.WALKS_OUT)
        self.assertEqual(result, expected)  # type: ignore

    def test_f_digest_index_gaponsetbase(self):
        """
        Tests GapOnSetBase is returned
        """
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCA-TAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCA-TAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]
        result = f_digest_index(self.create_array(seqs), self.config, 24, 0)

        # The Expected Sequence
        expected = (24, DIGESTION_ERROR.GAP_ON_SET_BASE)
        self.assertEqual(result, expected)  # type: ignore

    def test_f_digest_index_walktofar(self):
        """
        Tests WalksToFar is returned
        """

        seqs = [
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]

        local_cfg = Config()
        local_cfg.primer_max_walk = 10  # Force the primer to walk to far

        result = f_digest_index(self.create_array(seqs), local_cfg, 60, 0)

        # The Expected Sequence
        expected = (60, DIGESTION_ERROR.WALK_TO_FAR)
        self.assertEqual(result, expected)  # type: ignore


class TestDigest(unittest.TestCase):
    def setUp(self):
        self.config = Config()
        self.pm = PM()

    def create_array(self, seqs) -> np.ndarray:
        array_list = []
        for seq in seqs:
            array_list.append([x for x in seq])
        return np.array(array_list)

    def test_digestion_valid_fkmer(self):
        """Test the digestion"""
        self.maxDiff = None
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCANTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]
        msa_array = self.create_array(seqs)

        results = digest(
            msa_array=msa_array,
            config=self.config,
            indexes=([60], [1]),
            progress_manager=PM(),
        )
        expected_fkmer = FKmer([b"CCAATGGTGCAAAAGGTATAATCATTAAT"], 60)

        fkmer = results[0][0]
        self.assertEqual(fkmer.seqs(), expected_fkmer.seqs())
        self.assertEqual(fkmer.end, expected_fkmer.end)

    def test_digestion_valid_rkmer(self):
        """Test the digestion"""
        self.maxDiff = None
        seqs = [
            "CCAATGGTGCAAAAGGTATAATCANTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
            "CCAATGGTGCAAAAGGTATAATCATTAATGTCCAATGGTGCAAAAGGTATAATCATTAATGT",
        ]
        msa_array = msa_array = self.create_array(seqs)

        results = digest(
            msa_array=msa_array,
            config=self.config,
            indexes=([1], [25]),
            progress_manager=PM(),
        )
        expected_rkmer = RKmer([b"ATACCTTTTGCACCATTGGACATTA"], 25)

        rkmer = results[1][0]

        self.assertEqual(rkmer.seqs(), expected_rkmer.seqs())
        self.assertEqual(rkmer.start, expected_rkmer.start)


if __name__ == "__main__":
    unittest.main()
