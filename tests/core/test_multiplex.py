import pathlib
import unittest

import primalscheme3.core.config as config
from primalscheme3.core.classes import FKmer, MatchDB, PrimerPair, RKmer
from primalscheme3.core.multiplex import Multiplex


class TestMultiplex(unittest.TestCase):
    db_path = pathlib.Path("./tests/core/mulitplex").absolute()
    matchdb = MatchDB(db_path, [], 30)  # Create an empty matchdb

    # Create a config dict
    cfg = config.config_dict
    cfg["npools"] = 2
    cfg["mismatch_fuzzy"] = True
    cfg["mismatch_kmersize"] = 20
    cfg["mismatch_product_size"] = 200

    def test_next_pool_2_pool(self):
        """
        Test if calling the next_pool method returns the correct pool
        """
        self.cfg["npools"] = 2
        multiplex = Multiplex(cfg=self.cfg, matchDB=self.matchdb)
        current_pool = multiplex._current_pool
        next_pool = multiplex.next_pool()
        self.assertEqual(
            current_pool + 1,
            next_pool,
        )

    def test_next_pool_1_pool(self):
        """
        Test if calling the next_pool method returns the correct pool
        """
        self.cfg["npools"] = 1
        multiplex = Multiplex(cfg=self.cfg, matchDB=self.matchdb)
        current_pool = multiplex._current_pool
        next_pool = multiplex.next_pool()
        self.assertEqual(
            current_pool,
            next_pool,
        )

    def test_add_primer_pair_to_pool(self):
        """
        Test if method add_primer_pair_to_pool does whats expected
        """
        self.cfg["npools"] = 2
        multiplex = Multiplex(cfg=self.cfg, matchDB=self.matchdb)
        pp_msa_index = 0

        primerpair = PrimerPair(FKmer(10, "A"), RKmer(20, "T"), pp_msa_index)
        # Add a primerpair to pool 0
        multiplex.add_primer_pair_to_pool(primerpair, 0, pp_msa_index)

        # Check that the primerpair has beed added to _last_pp_added
        self.assertEqual(
            multiplex._last_pp_added[-1],
            primerpair,
        )
        # Check that the primerpair has been added to the correct pool
        self.assertEqual(multiplex._pools[0], [primerpair])
        # Check that current pool has updated
        self.assertEqual(multiplex._current_pool, 1)
        # Check that the primerpair has had its msa_index updated
        self.assertEqual(primerpair.msa_index, pp_msa_index)
        # Check amplicon number has been assigned
        self.assertEqual(multiplex._last_pp_added[-1].amplicon_number, 1)

    def test_remove_last_primer_pair(self):
        self.cfg["npools"] = 2
        self.cfg["minoverlap"] = 10
        multiplex = Multiplex(cfg=self.cfg, matchDB=self.matchdb)
        primerpair = PrimerPair(FKmer(100, "AA"), RKmer(200, "TT"), None)

        # Add a primerpair to pool 0
        multiplex.add_primer_pair_to_pool(primerpair, multiplex._current_pool, 0)

        # Remove the last primerpair
        last_pp = multiplex.remove_last_primer_pair()

        # Check that the lsat primerpair has been returned
        self.assertEqual(last_pp, primerpair)
        # Check the primer has been removed from the _last_pp_added
        self.assertEqual(len(multiplex._last_pp_added), 0)
        # Check the primer has been removed from the pool
        self.assertEqual(len(multiplex._pools[0]), 0)
        # Check the current pool has been reset
        self.assertEqual(multiplex._current_pool, 0)
        # Print last_pp has the expected pool
        self.assertEqual(last_pp.pool, 0)

    def test_does_overlap(self):
        """
        Test if method does_overlap does whats expected
        """
        self.cfg["npools"] = 2
        multiplex = Multiplex(cfg=self.cfg, matchDB=self.matchdb)

        # Create a primerpair
        primerpair = PrimerPair(FKmer(10, "A"), RKmer(200, "T"), 0)
        # Add a primerpair to pool 0
        multiplex.add_primer_pair_to_pool(primerpair, multiplex._current_pool, 0)

        # Check that the primerpair does overlap itself
        self.assertTrue(multiplex.does_overlap(primerpair, 0))

        # Check that nonoverlapping primerpair in the same msa do not does not overlap
        primerpair = PrimerPair(FKmer(1000, "A"), RKmer(2000, "T"), 0)
        self.assertFalse(multiplex.does_overlap(primerpair, 0))

        # Check that overlapping primerpair in differnt msas do not overlap
        primerpair_newmsa = PrimerPair(FKmer(100, "A"), RKmer(200, "T"), 1)
        self.assertFalse(multiplex.does_overlap(primerpair_newmsa, 0))

        # Check that an overlapping primerpair in a different pool does not overlap
        self.assertFalse(multiplex.does_overlap(primerpair, 1))

    def test_all_primerpairs(self):
        """
        Test if method all_primerpairs does whats expected
        """
        self.cfg["npools"] = 2
        multiplex = Multiplex(cfg=self.cfg, matchDB=self.matchdb)

        # Create a primerpair
        primerpair = PrimerPair(FKmer(10, "A"), RKmer(200, "T"), 0)
        # Add a primerpair to pool 0
        multiplex.add_primer_pair_to_pool(primerpair, multiplex._current_pool, 0)

        # Check that the primerpair does overlap itself
        self.assertEqual(multiplex.all_primerpairs(), [primerpair])


if __name__ == "__main__":
    unittest.main()
