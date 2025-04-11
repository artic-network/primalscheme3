import unittest

import numpy as np

from primalscheme3.core.mapping import (
    check_for_end_on_gap,
    create_mapping,
    fix_end_on_gap,
    generate_consensus,
    ref_index_to_msa,
)


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
            [
                None,
                None,
                None,
                0,
                1,
                2,
                3,
                4,
                None,
                None,
                5,
                6,
                7,
                8,
                9,
                None,
                None,
                None,
                None,
            ],
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


class Test_RefIndexToMsa(unittest.TestCase):
    def test_ref_index_to_msa(self):
        mapping_array = np.array(
            [
                None,
                None,
                None,
                0,
                1,
                2,
                3,
                4,
                None,
                None,
                5,
                6,
                7,
                8,
                9,
                None,
                None,
                None,
                None,
            ]
        )
        expected_answer = {
            0: 3,
            1: 4,
            2: 5,
            3: 6,
            4: 7,
            5: 10,
            6: 11,
            7: 12,
            8: 13,
            9: 14,
            10: 15,
        }
        result = ref_index_to_msa(mapping_array)
        self.assertEqual(result, expected_answer)


class Test_CheckForEndOnGap(unittest.TestCase):
    input = np.array(
        [
            [x for x in "---ATCGA--TCAGC----"],
            [x for x in "TTTATCGATTTCAGCACTG"],
            [x for x in "TTTATCGATTTCAGCACTG"],
            [x for x in "ATTATCGATTTCAGCACTG"],
        ]
    )
    mapping_array, array = create_mapping(input, 0)
    ref_index_to_msa_dict = ref_index_to_msa(mapping_array)

    def test_check_for_end_on_gap(self):
        result = check_for_end_on_gap(self.ref_index_to_msa_dict, 5)
        self.assertTrue(result)

    def test_check_for_end_on_gap_false(self):
        result = check_for_end_on_gap(self.ref_index_to_msa_dict, 3)
        self.assertFalse(result)


class Test_FixEndOnGap(unittest.TestCase):
    input = np.array(
        [
            [x for x in "---ATCGA--TCAGC----"],
            [x for x in "TTTATCGATTTCAGCACTG"],
            [x for x in "TTTATCGATTTCAGCACTG"],
            [x for x in "ATTATCGATTTCAGCACTG"],
        ]
    )
    mapping_array, array = create_mapping(input, 0)
    ref_index_to_msa_dict = ref_index_to_msa(mapping_array)

    def test_fix_end_on_gap(self):
        result = fix_end_on_gap(self.ref_index_to_msa_dict, 5)
        self.assertEqual(result, 8)


if __name__ == "__main__":
    unittest.main()
