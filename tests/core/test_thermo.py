import unittest

from primalscheme3.core.config import Config
from primalscheme3.core.thermo import THERMO_RESULT, calc_tm, gc, max_homo, thermo_check


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
            # Tm 59.9 °C with the digest's 0.8 M formamide
            "GGTATAGTCAATTGCCTAGTCTTCTGTTC": THERMO_RESULT.PASS,
            "GGGGGGGCCCCCCCC": THERMO_RESULT.HIGH_GC,
            "AAAATTTAATATATAT": THERMO_RESULT.LOW_GC,
            # Longer than before so its Tm clears the window with formamide
            "GTAACAGATACGTTGCAAAGTTTTTTTGACAG": THERMO_RESULT.MAX_HOMOPOLY,
            "AG": THERMO_RESULT.LOW_TM,
            "AGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA": THERMO_RESULT.TO_LONG,
        }

        for seq, truth in test_data.items():
            self.assertEqual(thermo_check(seq, config=self.config), truth, msg=seq)

    def test_tm_window_matches_the_digest(self):
        # Tm with the config's formamide, failing at or above the maximum (no
        # +2 °C slack), as the Rust digest does
        seq = "GGTATAGTCAATTGCCTAGTCTTCTGTTC"
        c = Config()
        tm = calc_tm(
            seq, c.mv_conc, c.dv_conc, c.dntp_conc, c.dna_conc, c.formamide_conc
        )
        self.assertEqual(
            thermo_check(seq, Config(primer_tm_max=tm + 0.01)), THERMO_RESULT.PASS
        )
        self.assertEqual(
            thermo_check(seq, Config(primer_tm_max=tm)), THERMO_RESULT.HIGH_TM
        )
        self.assertEqual(
            thermo_check(seq, Config(primer_tm_min=tm + 0.01)), THERMO_RESULT.LOW_TM
        )
        # Formamide lowers it, as in the digest
        self.assertLess(
            tm, calc_tm(seq, c.mv_conc, c.dv_conc, c.dntp_conc, c.dna_conc, 0.0) - 1.0
        )

    def test_wrappers_are_primalschemers(self):
        import primalschemers as ps

        from primalscheme3.core.thermo import calc_hairpin_tm, forms_hairpin

        c = Config()
        seq = "CCCCCATCCGATCAGGGGG"
        kw = dict(
            mv_conc=c.mv_conc,
            dv_conc=c.dv_conc,
            dntp_conc=c.dntp_conc,
            dna_conc=c.dna_conc,
        )
        hairpin = ps.calc_hairpin(seq, **kw).tm
        # primer3 2.3.0 calc_hairpin at these conditions (thal is identical)
        self.assertAlmostEqual(hairpin, 67.96550347919458, places=6)
        self.assertEqual(calc_hairpin_tm(seq, **kw), hairpin)
        self.assertTrue(forms_hairpin([seq], Config(primer_hairpin_th_max=60)))
        self.assertFalse(forms_hairpin([seq], Config(primer_hairpin_th_max=70)))


if __name__ == "__main__":
    unittest.main()
