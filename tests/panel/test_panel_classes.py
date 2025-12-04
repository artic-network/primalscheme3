import pathlib
import unittest

from primalscheme3.core.config import Config
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.panel.panel_classes import PanelMSA, does_overlap


class TestPanelMSA(unittest.TestCase):
    config = Config()

    def test_pointer(self):
        """
        Test case to check if the pointer is set correctly.
        """
        msa = PanelMSA(
            "test",
            pathlib.Path("tests/core/test_mismatch.fasta"),
            0,
            ProgressManager(),
            config=self.config,
        )
        # Add fake primerpairs to the MSA
        msa.primerpairs = [(x, 10, 0) for x in range(10)]  # type: ignore

        # Check if the pointer is set correctly
        self.assertEqual(msa.primerpairpointer, 0)

        # Assert that all primerpairs are returned
        self.assertEqual(msa.primerpairs, list(msa.iter_unchecked_primerpairs()))

        # Update the pointer
        msa.primerpairpointer = 5
        # Assert the correct primerpairs are returned
        self.assertEqual(list(msa.iter_unchecked_primerpairs()), msa.primerpairs[5:])


class TestDoesOverlap(unittest.TestCase):
    def test_no_overlap(self):
        """
        Test case to check if there is no overlap between new_pp and current_pps.
        """
        new_pp = (10, 20, 0)
        current_pps = [(0, 5, 0), (25, 30, 0)]
        self.assertFalse(does_overlap(new_pp, current_pps))

    def test_overlap(self):
        """
        Test case to check if there is overlap between new_pp and current_pps.
        """
        new_pp = (10, 20, 0)
        current_pps = [(0, 15, 0), (25, 30, 0)]
        self.assertTrue(does_overlap(new_pp, current_pps))

    def test_overlap_with_same_range(self):
        """
        Test case to check if there is overlap between new_pp and current_pps with the same range.
        """
        new_pp = (10, 20, 0)
        current_pps = [(10, 20, 0), (25, 30, 0)]
        self.assertTrue(does_overlap(new_pp, current_pps))

    def test_overlap_with_different_msa(self):
        """
        Test case to check if there is overlap between new_pp and current_pps with different msa.
        """
        new_pp = (10, 20, 0)
        current_pps = [(0, 15, 1), (25, 30, 1)]
        self.assertFalse(does_overlap(new_pp, current_pps))


if __name__ == "__main__":
    unittest.main()
