import pathlib
import unittest

import numpy as np
from Bio import SeqIO

from primalscheme3.core.errors import MSAFileInvalid, MSAFileInvalidLength
from primalscheme3.core.msa import parse_msa


class TestParseMSA(unittest.TestCase):
    def setUp(self):
        self.msa_diff_length = pathlib.Path(
            "tests/test_data/test_msa/test_msa_diff_length.fasta"
        )
        self.msa_id_collision = pathlib.Path(
            "tests/test_data/test_msa/test_msa_id_collision.fasta"
        )
        self.msa_non_fasta = pathlib.Path(
            "tests/test_data/test_msa/test_msa_non_fasta.fasta"
        )
        self.msa_empty_col = pathlib.Path(
            "tests/test_data/test_msa/test_msa_empty_col.fasta"
        )

    def test_parse_msa_diff_length(self):
        """
        Checks if the MSA contains sequences of different lengths
        """
        with self.assertRaises(MSAFileInvalidLength):
            _ = parse_msa(self.msa_diff_length)

    def test_parse_msa_id_collision(self):
        """
        Checks if the MSA contains two identical IDs
        """
        with self.assertRaises(MSAFileInvalid):
            _ = parse_msa(self.msa_id_collision)

    def test_parse_msa_non_fasta(self):
        """
        Checks if the MSA is not in fasta format
        """
        with self.assertRaises(MSAFileInvalid):
            _ = parse_msa(self.msa_non_fasta)

    def test_removes_empty_columns(self):
        """
        Checks if the MSA removes empty columns
        """
        # Get the original size of the arrays
        records_index = SeqIO.index(str(self.msa_empty_col), "fasta")
        array = np.array(
            [record.seq.upper() for record in records_index.values()],
            dtype="U1",
            ndmin=2,  # Enforce 2D array even if one genome
        )
        original_cols = array.shape[1]

        # parse the array
        new_array, _ = parse_msa(self.msa_empty_col)

        # Check if the new array is smaller than the original
        self.assertTrue(new_array.shape[1] == original_cols - 1)

        # Check the correct col was removed
        for col_index in range(0, new_array.shape[1]):
            slice: set[str] = set(new_array[:, col_index])
            self.assertTrue(slice != {"", "-"})


if __name__ == "__main__":
    unittest.main()
