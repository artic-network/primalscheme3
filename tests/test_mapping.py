import unittest
import numpy as np

from primalscheme3.core.mapping import trucnate_msa, create_mapping, generate_consensus


class Test_TruncateMsa(unittest.TestCase):
    input = np.array(
        [
            [
                "",
                "",
                "",
                "A",
                "T",
                "C",
                "T",
                "A",
                "-",
                "-",
                "T",
                "C",
                "A",
                "G",
                "C",
                "",
                "",
                "",
                "",
            ],
            [x for x in "TTTATCNATTTCAGCACTG"],
            [x for x in "TTTATCNATTTCAGCACTG"],
            [x for x in "ATTATCNATTTCAGCACTG"],
        ]
    )

    def test_trucnate_msa(self):
        """
        Test expand_ambs correctly expands ambiguity codes
        """
        test_input = self.input.copy()
        result = trucnate_msa(test_input, 0)

        # Check the result is as expected
        self.assertEqual(result.shape, (4, 12))


class Test_CreateMapping(unittest.TestCase):
    input = np.array(
        [
            [
                "",
                "",
                "",
                "A",
                "T",
                "C",
                "T",
                "A",
                "-",
                "-",
                "T",
                "C",
                "A",
                "G",
                "C",
                "",
                "",
                "",
                "",
            ],
            [x for x in "TTTATCNATTTCAGCACTG"],
            [x for x in "TTTATCNATTTCAGCACTG"],
            [x for x in "ATTATCNATTTCAGCACTG"],
        ]
    )

    def test_create_mapping(self):
        """
        Test expand_ambs correctly expands ambiguity codes
        """

        ## MAP  [0,    1,   2,   3,   4, None, None, 5,   6,   7,   8,   9]
        # trunc_array([
        #       ['A', 'T', 'C', 'T', 'A', '-', '-', 'T', 'C', 'A', 'G', 'C'],
        #       ['A', 'T', 'C', 'N', 'A', 'T', 'T', 'T', 'C', 'A', 'G', 'C'],
        #       ['A', 'T', 'C', 'N', 'A', 'T', 'T', 'T', 'C', 'A', 'G', 'C'],
        #       ['A', 'T', 'C', 'N', 'A', 'T', 'T', 'T', 'C', 'A', 'G', 'C']],
        #   dtype='<U1'))

        test_input = self.input.copy()
        mapping_array, trunc_msa = create_mapping(test_input, 0)

        # Check the result is as expected
        self.assertEqual(
            list(mapping_array),
            [0, 1, 2, 3, 4, None, None, 5, 6, 7, 8, 9],
        )


class Test_GenerateConsensus(unittest.TestCase):
    def test_generate_consensus(self):
        input = np.array(
            [
                [x for x in "---ATCGA--TCAGC----"],
                [x for x in "TTTATCGATTTCAGCACTG"],
                [x for x in "TTTATCGATTTCAGCACTG"],
                [x for x in "ATTATCGATTTCAGCACTG"],
            ]
        )
        expected_answer = "TTTATCGATTTCAGCACTG"
        result = generate_consensus(input)
        self.assertEqual(result, expected_answer)

    def test_generate_consensus_all_n(self):
        """Having N in all positions should enable N to appear in concensus"""
        input = np.array(
            [
                [x for x in "---ATCNA--TCAGC----"],
                [x for x in "TTTATCNATTTCAGCACTG"],
                [x for x in "TTTATCNATTTCAGCACTG"],
                [x for x in "ATTATCNATTTCAGCACTG"],
            ]
        )
        expected_answer = "TTTATCNATTTCAGCACTG"
        result = generate_consensus(input)
        self.assertEqual(result, expected_answer)

    def test_generate_consensus_not_all_n(self):
        """Having N in all but one positions should prevent N from appearing in concensus"""
        input = np.array(
            [
                [x for x in "---ATCTA--TCAGC----"],
                [x for x in "TTTATCNATTTCAGCACTG"],
                [x for x in "TTTATCNATTTCAGCACTG"],
                [x for x in "ATTATCNATTTCAGCACTG"],
            ]
        )
        expected_answer = "TTTATCTATTTCAGCACTG"
        result = generate_consensus(input)
        self.assertEqual(result, expected_answer)
