import unittest

from primalscheme3.interaction.interaction import (
    create_cigar,
    create_str,
)


class TestInteraction(unittest.TestCase):
    def test_create_cigar(self):
        """
        Test the create_cigar function with various sequence pairs
        """

        # Test complementary match
        seq1 = "ATGC"
        seq2 = "TACG"
        expected = "||||"
        self.assertEqual(create_cigar(seq1, seq2), expected)

        # Test partial match
        seq1 = "ATGC"
        seq2 = "TTGC"
        expected = "|..."
        self.assertEqual(create_cigar(seq1, seq2), expected)

        # Test with non-standard bases
        seq1 = "ATGCN"
        seq2 = "TACGN"
        expected = "|||| "
        self.assertEqual(create_cigar(seq1, seq2), expected)

    def test_create_str(self):
        """
        Test the create_str function with different offsets
        """
        seq1 = "ATGC"
        seq2 = "GCAT"

        # Test with no offset
        expected_no_offset = "score: 0.0\n5'-ATGC-3' >\n   ||||     \n3'-TACG-5'  \n"
        self.assertEqual(create_str(seq1, seq2, 0, 0.0), expected_no_offset)

        # Test with positive offset
        expected_pos_offset = (
            "score: 0.0\n 5'-ATGC-3' >\n    ...      \n3'-TACG-5'   \n"
        )
        self.assertEqual(create_str(seq1, seq2, 1, 0.0), expected_pos_offset)

        # Test with negative offset
        expected_neg_offset = "score: 0.0\n5'-ATGC-3' >\n    ...     \n 3'-TACG-5' \n"
        self.assertEqual(create_str(seq1, seq2, -1, 0.0), expected_neg_offset)


if __name__ == "__main__":
    unittest.main()
