import unittest

from primalscheme3.core.config import Config
from primalscheme3.core.thermo import THERMO_RESULT, gc, max_homo, thermo_check


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


class Test_PassesThermoChecks(unittest.TestCase):
    config = Config()

    def test_thermo_check(self):
        """
        Valuation order.
        """
        test_data = {
            "GTAATTCAGATACTGGTTGCAAAGTTATTATGA": THERMO_RESULT.PASS,
            "GGGGGGGCCCCCCCC": THERMO_RESULT.HIGH_GC,
            "AAAATTTAATATATAT": THERMO_RESULT.LOW_GC,
            "GTAACAGATACGTTGCAAAGTTTTTTTGA": THERMO_RESULT.MAX_HOMOPOLY,
            "AG": THERMO_RESULT.LOW_TM,
            "AGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA": THERMO_RESULT.TO_LONG,
        }

        for seq, truth in test_data.items():
            self.assertEqual(thermo_check(seq, config=self.config), truth, msg=seq)


if __name__ == "__main__":
    unittest.main()
