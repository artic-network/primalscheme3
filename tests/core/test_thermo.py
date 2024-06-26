import unittest

from primalscheme3.core.config import Config
from primalscheme3.core.thermo import THERMORESULT, gc, max_homo, passes_thermo_checks


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

    def test_passes_thermo_checks(self):
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
            self.assertEqual(passes_thermo_checks(seq, config=self.config), truth)


if __name__ == "__main__":
    unittest.main()
