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


if __name__ == "__main__":
    unittest.main()
