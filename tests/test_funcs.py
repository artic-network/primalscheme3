import unittest
import numpy as np

from primal_digest.seq_functions import (
    expand_ambs,
    remove_end_insertion,
    get_most_common_base,
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


if __name__ == "__main__":
    unittest.main()
