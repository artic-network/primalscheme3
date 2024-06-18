import pathlib
import unittest

import numpy as np

import primalscheme3.core.config as config
from primalscheme3.core.bedfiles import BedPrimerPair
from primalscheme3.core.classes import FKmer, MatchDB, PrimerPair, RKmer
from primalscheme3.core.msa import MSA
from primalscheme3.core.multiplex import Multiplex


class TestMultiplex(unittest.TestCase):
    db_path = pathlib.Path("./tests/core/mulitplex").absolute()
    matchdb = MatchDB(db_path, [], 30)  # Create an empty matchdb
    inputfile_path = pathlib.Path("./tests/core/test_mismatch.fasta").absolute()

    # Create a config dict
    cfg = config.config_dict
    cfg["npools"] = 2
    cfg["mismatch_fuzzy"] = True
    cfg["mismatch_kmersize"] = 20
    cfg["mismatch_product_size"] = 200

    # Create an MSA object
    msa = MSA("test", inputfile_path, 0, "first", None)

    def test_next_pool_2_pool(self):
        """
        Test if calling the next_pool method returns the correct pool
        """
        self.cfg["npools"] = 2
        multiplex = Multiplex(
            cfg=self.cfg, matchDB=self.matchdb, msa_dict={0: self.msa}
        )
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
        multiplex = Multiplex(
            cfg=self.cfg, matchDB=self.matchdb, msa_dict={0: self.msa}
        )
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
        multiplex = Multiplex(
            cfg=self.cfg, matchDB=self.matchdb, msa_dict={0: self.msa}
        )
        pp_msa_index = 0

        primerpair = PrimerPair(FKmer(10, ["A"]), RKmer(20, ["T"]), pp_msa_index)
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
        multiplex = Multiplex(
            cfg=self.cfg, matchDB=self.matchdb, msa_dict={0: self.msa}
        )
        primerpair = PrimerPair(FKmer(100, ["AA"]), RKmer(200, ["TT"]), 0)

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
        multiplex = Multiplex(
            cfg=self.cfg, matchDB=self.matchdb, msa_dict={0: self.msa, 1: self.msa}
        )

        # Create a primerpair
        primerpair = PrimerPair(FKmer(10, ["A"]), RKmer(200, ["T"]), 0)
        # Add a primerpair to pool 0
        multiplex.add_primer_pair_to_pool(primerpair, multiplex._current_pool, 0)

        # Check that the primerpair does overlap itself
        self.assertTrue(multiplex.does_overlap(primerpair, 0))

        # Check that nonoverlapping primerpair in the same msa do not does not overlap
        primerpair = PrimerPair(FKmer(1000, ["A"]), RKmer(2000, ["T"]), 0)
        self.assertFalse(multiplex.does_overlap(primerpair, 0))

        # Check that overlapping primerpair in differnt msas do not overlap
        primerpair_newmsa = PrimerPair(FKmer(100, ["A"]), RKmer(200, ["T"]), 1)
        self.assertFalse(multiplex.does_overlap(primerpair_newmsa, 0))

        # Check that an overlapping primerpair in a different pool does not overlap
        self.assertFalse(multiplex.does_overlap(primerpair, 1))

    def test_all_primerpairs(self):
        """
        Test if method all_primerpairs does whats expected
        """
        self.cfg["npools"] = 2
        multiplex = Multiplex(
            cfg=self.cfg, matchDB=self.matchdb, msa_dict={0: self.msa}
        )

        # Create a primerpair
        primerpair = PrimerPair(FKmer(10, ["A"]), RKmer(200, ["T"]), 0)
        # Add a primerpair to pool 0
        multiplex.add_primer_pair_to_pool(primerpair, multiplex._current_pool, 0)

        # Check that the primerpair does overlap itself
        self.assertEqual(multiplex.all_primerpairs(), [primerpair])

    def test_coverage(self):
        """
        Test if method coverage does whats expected
        """

        self.cfg["npools"] = 2
        multiplex = Multiplex(
            cfg=self.cfg, matchDB=self.matchdb, msa_dict={0: self.msa}
        )

        # Create a primerpair
        primerpair = PrimerPair(FKmer(10, ["A"]), RKmer(200, ["T"]), 0)

        # Check coverage is created correctly
        self.assertEqual(multiplex._coverage[0].sum(), 0)

        # Add a primerpair
        multiplex.update_coverage(0, primerpair, add=True)

        # Check that the primerpair coverage has been added
        self.assertEqual(multiplex._coverage[0].sum(), 190)

        # Remove the primerpair
        multiplex.update_coverage(0, primerpair, add=False)

        # Check that the primerpair coverage has been removed
        self.assertEqual(multiplex._coverage[0].sum(), 0)

    def test_coverage_circular(self):
        """
        Test if method coverage does whats expected
        """

        self.cfg["npools"] = 2
        multiplex = Multiplex(
            cfg=self.cfg, matchDB=self.matchdb, msa_dict={0: self.msa}
        )

        # Create a primerpair
        primerpair = PrimerPair(
            FKmer(len(self.msa.array[0]) - 100, ["A"]), RKmer(10, ["T"]), 0
        )

        # Check coverage is created correctly
        self.assertEqual(multiplex._coverage[0].sum(), 0)

        # Add a primerpair
        multiplex.update_coverage(0, primerpair, add=True)

        # Check that the primerpair coverage has been added
        self.assertEqual(multiplex._coverage[0].sum(), 110)

        # Remove the primerpair
        multiplex.update_coverage(0, primerpair, add=False)

        # Check that the primerpair coverage has been removed
        self.assertEqual(multiplex._coverage[0].sum(), 0)

    def test_lookup(self):
        """
        Test if method lookup does whats expected
        """

        self.cfg["npools"] = 2
        multiplex = Multiplex(
            cfg=self.cfg, matchDB=self.matchdb, msa_dict={0: self.msa}
        )

        msa_index = 0
        pool = 0
        # Create a primerpair
        primerpair = PrimerPair(FKmer(10, ["A"]), RKmer(200, ["T"]), 0)
        primerpair.pool = pool

        # Check lookup is created empty, in the correct shape
        for p in range(0, self.cfg["npools"]):
            self.assertEqual(np.count_nonzero(multiplex._lookup[msa_index][p, :]), 0)

        # Add a primerpair
        multiplex.update_lookup(primerpair, add=True)

        # Check that the primerpair coverage has been added
        self.assertEqual(np.count_nonzero(multiplex._lookup[msa_index][pool, :]), 192)

        # Remove the primerpair
        multiplex.update_lookup(primerpair, add=False)

        # Check that the primerpair coverage has been removed
        self.assertEqual(
            np.count_nonzero(multiplex._lookup[msa_index]),
            0,
        )

    def test_lookup_circular(self):
        """
        Test if method lookup does whats expected
        """

        self.cfg["npools"] = 2
        multiplex = Multiplex(
            cfg=self.cfg, matchDB=self.matchdb, msa_dict={0: self.msa}
        )

        msa_index = 0
        pool = 0
        # Create a primerpair
        primerpair = PrimerPair(
            FKmer(len(self.msa.array[0]) - 100, ["A"]), RKmer(10, ["T"]), 0
        )
        primerpair.pool = pool

        # Check lookup is created empty, in the correct shape
        for p in range(0, self.cfg["npools"]):
            self.assertEqual(np.count_nonzero(multiplex._lookup[msa_index][p, :]), 0)

        # Add a primerpair
        multiplex.update_lookup(primerpair, add=True)

        # Check that the primerpair coverage has been added
        self.assertEqual(np.count_nonzero(multiplex._lookup[msa_index][pool, :]), 112)

        # Remove the primerpair
        multiplex.update_lookup(primerpair, add=False)

        # Check that the primerpair coverage has been removed
        self.assertEqual(
            np.count_nonzero(multiplex._lookup[msa_index]),
            0,
        )

    def test_bedprimer(self):
        self.cfg["npools"] = 2
        multiplex = Multiplex(
            cfg=self.cfg, matchDB=self.matchdb, msa_dict={0: self.msa}
        )
        # Create a primerpair
        bedprimerpair = BedPrimerPair(
            FKmer(10, ["A"]),
            RKmer(100, ["T"]),
            msa_index=-1,
            chrom_name="test",
            amplicon_prefix="test",
            amplicon_number=1,
            pool=0,
        )

        # Add a primerpair to pool 0
        multiplex.add_primer_pair_to_pool(bedprimerpair, multiplex._current_pool, 0)

        # Check that the primerpair has beed added to _last_pp_added
        self.assertEqual(
            multiplex._last_pp_added[-1],
            bedprimerpair,
        )

        # Check that the primerpair has been added to the correct pool
        self.assertEqual(multiplex._pools[0], [bedprimerpair])

        # check the coverage has not changed
        self.assertEqual(multiplex._coverage[0].sum(), 0)
        # Check the lookup has not changed
        self.assertEqual(
            np.count_nonzero(multiplex._lookup[0][0, :]),
            0,
        )

        # Check primerpair can be removed
        multiplex.remove_last_primer_pair()
        self.assertEqual(len(multiplex._last_pp_added), 0)


if __name__ == "__main__":
    unittest.main()
