import pathlib
import unittest

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
        with self.assertRaises(MSAFileInvalidLength):
            _ = parse_msa(self.msa_id_collision)

    def test_parse_msa_non_fasta(self):
        """
        Checks if the MSA is not in fasta format
        """
        with self.assertRaises(MSAFileInvalid):
            _ = parse_msa(self.msa_non_fasta)


if __name__ == "__main__":
    unittest.main()
