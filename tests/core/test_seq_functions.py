import unittest

import numpy as np

from primalscheme3.core.seq_functions import (
    calc_entropy,
    calc_probs,
    expand_ambs,
    extend_ambiguous_base,
    get_most_common_base,
    remove_end_insertion,
    reverse_complement,
)


class Test_ExpandAmbs(unittest.TestCase):
    def test_expand_ambs(self):
        """
        Test expand_ambs correctly expands ambiguity codes
        """
        seqeuence = {"ATGM"}
        result = expand_ambs(seqeuence)
        self.assertEqual(result, {"ATGC", "ATGA"})

    def test_expand_ambs_multi(self):
        """
        Test expand_ambs correctly expands ambiguity codes, on mutliple seqs
        """
        seqeuence = {"ATGM", "ATGB"}
        result = expand_ambs(seqeuence)
        self.assertEqual(result, {"ATGC", "ATGA", "ATGT", "ATGG"})


class Test_RemoveEndInsertion(unittest.TestCase):
    def test_remove_end_insertion(self):
        """
        Ensure remove_end_insertion removes whats expected
        """
        input = np.array(
            [[x for x in "---ATCGA--TCAGC----"], [x for x in "TTTATCGATTTCAGCACTG"]]
        )
        expected_answer = {"ATCGA--TCAGC", "TTTATCGATTTCAGCACTG"}

        result = remove_end_insertion(input)
        self.assertEqual({"".join(x) for x in result}, expected_answer)

    def test_remove_end_insertion_no_change(self):
        """
        Ensure remove_end_insertion removes whats expected
        """
        input = np.array(
            [[x for x in "TTTATCGATCAGCACTG"], [x for x in "TTTATCGATCAGCACTG"]]
        )
        expected_answer = {"TTTATCGATCAGCACTG"}

        result = remove_end_insertion(input)
        self.assertEqual({"".join(x) for x in result}, expected_answer)


class Test_GetMostCommonBase(unittest.TestCase):
    def test_get_most_common_base(self):
        input = np.array(
            [
                [x for x in "---ATCGA--TCAGC----"],
                [x for x in "TTTATCGATTTCAGCACTG"],
                [x for x in "TTTATCGATTTCAGCACTG"],
                [x for x in "ATTATCGATTTCAGCACTG"],
            ]
        )
        expected_answer = "T"
        result = get_most_common_base(input, 0)
        self.assertEqual(result, expected_answer)


class Test_ReverseComplement(unittest.TestCase):
    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("ATCG"), "CGAT")
        self.assertEqual(reverse_complement("atcg"), "CGAT")
        self.assertEqual(reverse_complement("AGCT"), "AGCT")
        self.assertEqual(reverse_complement("NNNN"), "NNNN")
        self.assertEqual(reverse_complement("ACGTN"), "NACGT")

    def test_reverse_complement_error(self):
        """
        Should raise keys error if invalid base is passed
        """
        self.assertRaises(KeyError, reverse_complement, ".")


class Test_ExtendAmbiguousBase(unittest.TestCase):
    def test_extend_ambiguous_base(self):
        self.assertEqual(extend_ambiguous_base("M"), ["A", "C"])
        self.assertEqual(extend_ambiguous_base("R"), ["A", "G"])
        self.assertEqual(extend_ambiguous_base("W"), ["A", "T"])
        self.assertEqual(extend_ambiguous_base("S"), ["C", "G"])
        self.assertEqual(extend_ambiguous_base("Y"), ["C", "T"])
        self.assertEqual(extend_ambiguous_base("K"), ["G", "T"])
        self.assertEqual(extend_ambiguous_base("V"), ["A", "C", "G"])
        self.assertEqual(extend_ambiguous_base("H"), ["A", "C", "T"])
        self.assertEqual(extend_ambiguous_base("D"), ["A", "G", "T"])
        self.assertEqual(extend_ambiguous_base("B"), ["C", "G", "T"])

    def test_extend_ambiguous_base_error(self):
        """
        Invalid base should return "N"
        """
        self.assertEqual(extend_ambiguous_base("."), ["N"])


class TestCalcEntropy(unittest.TestCase):
    def test_calc_entropy(self):
        """
        Test calc_entropy returns the correct entropy value
        """
        probs = [0.5, 0.25, 0.125, 0.125]
        result = calc_entropy(probs)
        self.assertAlmostEqual(result, 1.75)

    def test_calc_entropy_zero_prob(self):
        """
        Test calc_entropy handles zero probability values
        """
        probs = [0.5, 0.25, 0.125, 0.0, 0.125]
        result = calc_entropy(probs)
        self.assertAlmostEqual(result, 1.75)

    def test_calc_entropy_all_zero_prob(self):
        """
        Test calc_entropy handles all zero probability values
        """
        probs = [0.0, 0.0, 0.0, 0.0]
        result = calc_entropy(probs)
        self.assertAlmostEqual(result, 0.0)

    def test_calc_entropy_single_prob(self):
        """
        Test calc_entropy handles single probability value
        """
        probs = [1.0]
        result = calc_entropy(probs)
        self.assertAlmostEqual(result, 0.0)


class TestCalcProbs(unittest.TestCase):
    def test_calc_probs_base(self):
        """
        Test calc_probs handles ambiguous bases correctly
        """
        bases = ["A", "A", "A", "A"]

        expected_result = np.array([1, 1, 1, 1])
        result = np.array(calc_probs(bases))

        self.assertTrue(np.isclose(result, expected_result).all())


if __name__ == "__main__":
    unittest.main()
