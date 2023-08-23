import unittest
from primal_digest.get_window import (
    get_f_window_FAST2,
    get_r_window_FAST2,
    get_pp_window,
)
from primal_digest.classes import FKmer, RKmer, PrimerPair


class Test_GetFWindowFAST2(unittest.TestCase):
    def test_get_f_window_FAST2(self):
        fkmers = [FKmer(end, "AA") for end in range(20, 100)]

        # Get all kmers that start between 40 and 50
        expected_fkmers = [
            fkmer for fkmer in fkmers if fkmer.end >= 40 and fkmer.end <= 50
        ]
        result_fkmers = get_f_window_FAST2(fkmers, 40, 50)
        self.assertEqual(result_fkmers, expected_fkmers)


class Test_GetRWindowFAST2(unittest.TestCase):
    def test_get_f_window_FAST2(self):
        rkmers = [RKmer(start, "AA") for start in range(20, 100)]

        # Get all kmers that start between 40 and 50
        expected_rkmers = [
            rkmer for rkmer in rkmers if rkmer.start >= 40 and rkmer.start <= 50
        ]
        result_rkmers = get_r_window_FAST2(rkmers, 40, 50)
        self.assertEqual(result_rkmers, expected_rkmers)


class Test_GetPpWindow(unittest.TestCase):
    def test_get_pp_window_ol(self):
        fkmers = [FKmer(end, "A") for end in range(20, 30)]
        rkmers = [RKmer(start, "A") for start in range(90, 100)]

        primerpairs = []
        for fkmer in fkmers:
            for rkmer in rkmers:
                primerpairs.append(PrimerPair(fkmer, rkmer, 0))
        primerpairs.sort(key=lambda pp: (pp.fprimer.end, pp.rprimer.start))

        # Get all kmers that start between 40 and 50
        expected_pos_ol_pp = [
            pp
            for pp in primerpairs
            if pp.fprimer.end >= 22 and pp.fprimer.end <= 28 and pp.rprimer.start >= 95
        ]

        result_pp = get_pp_window(primerpairs, 22, 28, 95)
        self.assertEqual(expected_pos_ol_pp, result_pp)


if __name__ == "__main__":
    unittest.main()
