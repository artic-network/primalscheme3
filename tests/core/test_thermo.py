import unittest

from primalscheme3.core.config import Config
from primalscheme3.core.thermo import (
    THERMORESULT,
    forms_hairpin,
    gc,
    max_homo,
    thermo_check,
)


class Test_GC(unittest.TestCase):
    def test_gc(self):
        """
        Test gc correctly calculates gc content
        """
        test_data = {
            "ACTGACTGC": 55.6,
            "GCTAGCTAGCTAGCTAGCTGATCGATCGT": 51.7,
            "GGGGGGGGGGGGGG": 100,
        }

        for seq, gc_truth in test_data.items():
            self.assertEqual(gc(seq), gc_truth)


class Test_MaxHomo(unittest.TestCase):
    def test_max_homo(self):
        test_data = {"ACGATCGATCGTAGCTTATCGAC": 2, "AAAA": 4, "ATCGTTTTTTTTTT": 10}

        for seq, mh_truth in test_data.items():
            self.assertEqual(max_homo(seq), mh_truth)


class Test_PassesThermoCHecks(unittest.TestCase):
    config = Config()

    def test_thermo_check(self):
        """
        Valution order.
        """
        test_data = {
            "GTAATTCAGATACTGGTTGCAAAGTTTTTATGA": THERMORESULT.PASS,
            "GGGGGGGCCCCCCCC": THERMORESULT.HIGH_GC,
            "AAAATTTAATATATAT": THERMORESULT.LOW_GC,
            "GTAATTCAGATACTGGTTGCAAAGTTTTTTTGA": THERMORESULT.MAX_HOMOPOLY,
            "AG": THERMORESULT.LOW_TM,
            "AGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA": THERMORESULT.HIGH_TM,
        }

        for seq, truth in test_data.items():
            self.assertEqual(thermo_check(seq, config=self.config), truth)


class Test_Forms_Hairpin(unittest.TestCase):
    config = Config()

    def test_no_hairpin(self):
        # Non hairpin seq
        seq = "CTCTTGTAGATCTGTTCTCTAAACGAACTTT"
        self.assertFalse(forms_hairpin(seq, self.config))

    def test_hairpin_mismatches(self):
        # Hairpin seq with 3' match
        seq = "GGGGGGGTAGATCTGTTCTCTAAACGCCCCC"
        self.assertTrue(forms_hairpin(seq, self.config))

        # Test that adding a two 3' mismatch will prevent hairpin formation
        mismatch_single_3p = seq + "T"
        self.assertFalse(forms_hairpin(mismatch_single_3p, self.config))

        # Test that adding a two 3' mismatch will prevent hairpin formation
        mismatch_double_3p = seq + "TT"
        self.assertFalse(forms_hairpin(mismatch_double_3p, self.config))

        # Test a 5' mismatch still triggers hairpin
        mismatch_5p = "T" + seq
        self.assertTrue(forms_hairpin(mismatch_5p, self.config))


if __name__ == "__main__":
    unittest.main()
