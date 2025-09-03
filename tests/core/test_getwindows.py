import unittest

from primalschemers import FKmer, RKmer  # type: ignore

from primalscheme3.core.classes import PrimerPair
from primalscheme3.core.get_window import (
    get_f_window_FAST2,
    get_pp_window,
    get_r_window_FAST2,
)


class Test_GetFWindowFAST2(unittest.TestCase):
    def test_get_f_window_FAST2(self):
        fkmers = [FKmer([b"AA"], end) for end in range(20, 100)]

        # Get all kmers that start between 40 and 50
        expected_fkmers = [
            fkmer for fkmer in fkmers if fkmer.end >= 40 and fkmer.end <= 50
        ]
        result_fkmers = get_f_window_FAST2(fkmers, 40, 50)
        self.assertEqual(result_fkmers, expected_fkmers)


class Test_GetRWindowFAST2(unittest.TestCase):
    def test_get_f_window_FAST2(self):
        rkmers = [RKmer([b"AA"], start) for start in range(20, 100)]

        # Get all kmers that start between 40 and 50
        expected_rkmers = [
            rkmer for rkmer in rkmers if rkmer.start >= 40 and rkmer.start <= 50
        ]
        result_rkmers = get_r_window_FAST2(rkmers, 40, 50)
        self.assertEqual(result_rkmers, expected_rkmers)


class Test_GetPpWindow(unittest.TestCase):
    def test_get_pp_window_ol(self):
        fkmers = [FKmer([b"A"], end) for end in range(20, 110)]
        rkmers = [RKmer([b"A"], start) for start in range(10, 100)]

        cfg = {"amplicon_size_min": 10, "amplicon_size_max": 70}
        msa_index = 0

        ## Generate all the primerpairs
        non_checked_pp = []
        for fkmer in fkmers:
            fkmer_start = min(fkmer.starts())
            # Get all rkmers that would make a valid amplicon
            pos_rkmer = get_r_window_FAST2(
                kmers=rkmers,
                start=fkmer_start + cfg["amplicon_size_min"],
                end=fkmer_start + cfg["amplicon_size_max"],
            )
            for rkmer in pos_rkmer:
                non_checked_pp.append(PrimerPair(fkmer, rkmer, msa_index))
        non_checked_pp.sort(key=lambda pp: (pp.fprimer.end, -pp.rprimer.start))

        fp_end_min = 22
        fp_end_max = 28
        rp_start_min = 95

        # Get all kmers that start between 40 and 50
        expected_pos_ol_pp = [
            pp
            for pp in non_checked_pp
            if pp.fprimer.end >= fp_end_min
            and pp.fprimer.end <= fp_end_max
            and pp.rprimer.start >= rp_start_min
        ]

        result_pp = get_pp_window(non_checked_pp, fp_end_min, fp_end_max, rp_start_min)
        self.assertEqual(expected_pos_ol_pp, result_pp)


if __name__ == "__main__":
    unittest.main()
