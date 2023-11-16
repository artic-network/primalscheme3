import unittest

from primalscheme3.panel.minimal_scheme_classes import does_overlap


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
