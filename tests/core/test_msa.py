import pathlib
import unittest

from primalscheme3.core.config import Config
from primalscheme3.core.errors import (
    MSAFileInvalid,
    MSAFileInvalidBase,
    MSAFileInvalidLength,
)
from primalscheme3.core.msa import MSA, parse_msa
from primalscheme3.core.progress_tracker import ProgressManager


class TestMSA(unittest.TestCase):
    msa_path = pathlib.Path("tests/test_data/test_msa/test_msa_valid.fasta").absolute()
    config = Config()

    def test_create_msa(self):
        # Should not raise anything
        pm = ProgressManager()
        _ = MSA(
            name="test",
            path=self.msa_path,
            msa_index=0,
            progress_manager=pm,
            config=self.config,
        )


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
        self.msa_non_dna = pathlib.Path(
            "tests/test_data/test_msa/test_msa_non_dna.fasta"
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

    def test_parse_msa_non_dna(self):
        """
        Checks if the MSA contains non-DNA characters
        """
        with self.assertRaises(MSAFileInvalidBase):
            _ = parse_msa(self.msa_non_dna)

    def test_valid_msa(self):
        """
        Checks if the MSA is valid
        """
        _ = parse_msa(pathlib.Path("tests/test_data/test_msa/test_msa_valid.fasta"))


if __name__ == "__main__":
    unittest.main()
