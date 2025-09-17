import pathlib
import unittest

import numpy as np
from primalschemers import FKmer, RKmer  # type: ignore

from primalscheme3.core.bedfiles import BedPrimerPair
from primalscheme3.core.classes import PrimerPair
from primalscheme3.core.config import Config
from primalscheme3.core.mismatches import MatchDB
from primalscheme3.core.msa import MSA
from primalscheme3.core.multiplex import Multiplex, PrimerPairCheck


class TestMultiplex(unittest.TestCase):
    db_path = pathlib.Path("./tests/core/mulitplex").absolute()
    config = Config()
    matchdb = MatchDB(db_path, [], config)  # Create an empty matchdb
    inputfile_path = pathlib.Path("./tests/core/test_mismatch.fasta").absolute()
    nCoV_2019_76_RIGHT_0 = "ACACCTGTGCCTGTTAAACCAT"
    nCoV_2019_18_LEFT_0 = "TGGAAATACCCACAAGTTAATGGTTTAAC"

    # Create a config dict
    config = Config()
    config.n_pools = 2

    # Create an MSA object
    msa = MSA("test", inputfile_path, 0, None, config)

    def test_next_pool_2_pool(self):
        """
        Test if calling the next_pool method returns the correct pool
        """
        self.config.n_pools = 2
        multiplex = Multiplex(
            config=self.config, matchDB=self.matchdb, msa_dict={0: self.msa}
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
        self.config.n_pools = 1
        multiplex = Multiplex(
            config=self.config, matchDB=self.matchdb, msa_dict={0: self.msa}
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
        self.config.n_pools = 2
        multiplex = Multiplex(
            config=self.config, matchDB=self.matchdb, msa_dict={0: self.msa}
        )
        pp_msa_index = 0

        primerpair = PrimerPair(FKmer([b"A"], 10), RKmer([b"T"], 20), pp_msa_index)
        # Add a primerpair to pool 0
        multiplex.add_primer_pair_to_pool(primerpair, 0, pp_msa_index)

        # Check that the primerpair has been added to _last_pp_added
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

    def test_add_primer_pair_to_pool_set_amp_number(self):
        """
        Test if method add_primer_pair_to_pool does what expected, with a specified amp number
        """
        pool = 1
        amp_num = 10

        self.config.n_pools = 2
        multiplex = Multiplex(
            config=self.config, matchDB=self.matchdb, msa_dict={0: self.msa}
        )
        pp_msa_index = 0

        primerpair = PrimerPair(FKmer([b"A"], 10), RKmer([b"T"], 20), pp_msa_index)
        primerpair.amplicon_number = amp_num
        # Add a primerpair to pool 1
        multiplex.add_primer_pair_to_pool(primerpair, pool, pp_msa_index)

        # Check that the primerpair has been added to _last_pp_added
        self.assertEqual(
            multiplex._last_pp_added[-1],
            primerpair,
        )
        # Check that the primerpair has been added to the correct pool
        self.assertEqual(multiplex._pools[pool], [primerpair])
        # Check that current pool has updated
        self.assertEqual(multiplex._current_pool, 0)
        # Check that the primerpair has had its msa_index updated
        self.assertEqual(primerpair.msa_index, pp_msa_index)
        # Check amplicon number has not changed
        self.assertEqual(multiplex._last_pp_added[-1].amplicon_number, amp_num)

    def test_remove_last_primer_pair(self):
        self.config.n_pools = 2
        multiplex = Multiplex(
            config=self.config, matchDB=self.matchdb, msa_dict={0: self.msa}
        )
        primerpair = PrimerPair(FKmer([b"AA"], 100), RKmer([b"TT"], 200), 0)

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
        self.config.n_pools = 2
        multiplex = Multiplex(
            config=self.config,
            matchDB=self.matchdb,
            msa_dict={0: self.msa, 1: self.msa},
        )

        # Create a primerpair
        primerpair = PrimerPair(FKmer([b"A"], 100), RKmer([b"T"], 200), 0)
        # Add a primerpair to pool 0
        multiplex.add_primer_pair_to_pool(primerpair, multiplex._current_pool, 0)

        # Check that the primerpair does overlap itself
        self.assertTrue(multiplex.does_overlap(primerpair, 0))

        # Check that very non-overlapping primerpair in the same msa do not does not overlap
        self.assertFalse(
            multiplex.does_overlap(
                PrimerPair(FKmer([b"A"], 300), RKmer([b"T"], 400), 0), 0
            )
        )

        # Check that overlapping primerpair in differnt msas do not overlap
        self.assertFalse(
            multiplex.does_overlap(
                PrimerPair(FKmer([b"A"], 100), RKmer([b"T"], 200), 1), 0
            )
        )

        # Check that an overlapping primerpair in a different pool does not overlap
        self.assertFalse(multiplex.does_overlap(primerpair, 1))

        # Check that circular overlapping primerpair overlap
        self.assertTrue(
            multiplex.does_overlap(
                PrimerPair(FKmer([b"A"], 900), RKmer([b"T"], 200), 0), 0
            )
        )

        # Check that circular non-overlapping primerpair doesn't overlap
        primerpair_circular = PrimerPair(FKmer([b"A"], 900), RKmer([b"T"], 50), 0)
        self.assertFalse(multiplex.does_overlap(primerpair_circular, 0))

        # Check for off-by-one error with linear overlapping primerpair
        #      100  200
        #       A    T
        self.assertFalse(
            multiplex.does_overlap(
                PrimerPair(FKmer([b"A"], 202), RKmer([b"T"], 300), 0), 0
            )
        )

    def test_all_primerpairs(self):
        """
        Test if method all_primerpairs does whats expected
        """
        self.config.n_pools = 2
        multiplex = Multiplex(
            config=self.config, matchDB=self.matchdb, msa_dict={0: self.msa}
        )

        # Create a primerpair
        primerpair = PrimerPair(FKmer([b"A"], 10), RKmer([b"T"], 200), 0)
        # Add a primerpair to pool 0
        multiplex.add_primer_pair_to_pool(primerpair, multiplex._current_pool, 0)

        # Check that the primerpair does overlap itself
        self.assertEqual(multiplex.all_primerpairs(), [primerpair])

    def test_coverage(self):
        """
        Test if method coverage does whats expected
        """

        self.config.n_pools = 2
        multiplex = Multiplex(
            config=self.config, matchDB=self.matchdb, msa_dict={0: self.msa}
        )

        # Create a primerpair
        primerpair = PrimerPair(FKmer([b"A"], 10), RKmer([b"T"], 200), 0)

        # Check coverage is created correctly
        self.assertEqual(multiplex._coverage[0].sum(), 0)

        # Add a primerpair
        multiplex.update_coverage(primerpair, add=True)

        # Check that the primerpair coverage has been added
        self.assertEqual(multiplex._coverage[0].sum(), 190)

        # Remove the primerpair
        multiplex.update_coverage(primerpair, add=False)

        # Check that the primerpair coverage has been removed
        self.assertEqual(multiplex._coverage[0].sum(), 0)

    def test_coverage_circular(self):
        """
        Test if method coverage does whats expected
        """

        self.config.n_pools = 2
        multiplex = Multiplex(
            config=self.config, matchDB=self.matchdb, msa_dict={0: self.msa}
        )

        # Create a primerpair
        primerpair = PrimerPair(
            FKmer([b"A"], len(self.msa.array[0]) - 100), RKmer([b"T"], 10), 0
        )

        # Check coverage is created correctly
        self.assertEqual(multiplex._coverage[0].sum(), 0)

        # Add a primerpair
        multiplex.update_coverage(primerpair, add=True)

        # Check that the primerpair coverage has been added
        self.assertEqual(multiplex._coverage[0].sum(), 110)

        # Remove the primerpair
        multiplex.update_coverage(primerpair, add=False)

        # Check that the primerpair coverage has been removed
        self.assertEqual(multiplex._coverage[0].sum(), 0)

    def test_lookup(self):
        """
        Test if method lookup does whats expected
        """

        self.config.n_pools = 2
        multiplex = Multiplex(
            config=self.config, matchDB=self.matchdb, msa_dict={0: self.msa}
        )

        msa_index = 0
        pool = 0
        # Create a primerpair
        primerpair = PrimerPair(FKmer([b"A"], 10), RKmer([b"T"], 200), 0)
        primerpair.pool = pool

        # Check lookup is created empty, in the correct shape
        for p in range(0, self.config.n_pools):
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

        self.config.n_pools = 2
        multiplex = Multiplex(
            config=self.config, matchDB=self.matchdb, msa_dict={0: self.msa}
        )

        msa_index = 0
        pool = 0
        # Create a primerpair
        primerpair = PrimerPair(
            FKmer([b"A"], len(self.msa.array[0]) - 100), RKmer([b"T"], 10), 0
        )
        primerpair.pool = pool

        # Check lookup is created empty, in the correct shape
        for p in range(0, self.config.n_pools):
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
        self.config.n_pools = 2
        multiplex = Multiplex(
            config=self.config, matchDB=self.matchdb, msa_dict={0: self.msa}
        )
        # Create a primerpair
        bedprimerpair = BedPrimerPair(
            FKmer([b"A"], 10),
            RKmer([b"T"], 100),
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

    def test_check_add_primer(self):
        self.config.n_pools = 2
        multiplex = Multiplex(
            config=self.config,
            matchDB=self.matchdb,
            msa_dict={0: self.msa, 1: self.msa},
        )

        msa_index = 0
        pool = 0
        # Create a primerpair
        primerpair = PrimerPair(
            FKmer([self.nCoV_2019_18_LEFT_0.encode()], 50), RKmer([b"TA"], 100), 0
        )
        primerpair.pool = pool

        # Add the primerpair
        multiplex.add_primer_pair_to_pool(primerpair, msa_index, pool)

        # Check the same primerpair cannot be added again  due to overlap
        self.assertEqual(
            multiplex.check_primerpair_can_be_added(primerpair, pool),
            PrimerPairCheck.OVERLAP,
        )

        # Check a primerpair with a different msa_index can be added
        primerpair2 = PrimerPair(FKmer([b"TA"], 50), RKmer([b"TA"], 100), 1)
        self.assertEqual(
            multiplex.check_primerpair_can_be_added(primerpair2, pool),
            PrimerPairCheck.OK,
        )

        # Check a primerpair with a different pool can be added
        self.assertEqual(
            multiplex.check_primerpair_can_be_added(primerpair, pool + 1),
            PrimerPairCheck.OK,
        )

        # Check an interacting primerpair cannot be added
        interacting_primerpair = PrimerPair(
            FKmer([self.nCoV_2019_76_RIGHT_0.encode()], 150), RKmer([b"TA"], 200), 0
        )
        self.assertEqual(
            multiplex.check_primerpair_can_be_added(interacting_primerpair, pool),
            PrimerPairCheck.INTERACTING,
        )

    def test_get_next_amplicon_number(self):
        multiplex = Multiplex(
            config=self.config,
            matchDB=self.matchdb,
            msa_dict={0: self.msa, 1: self.msa},
        )
        msa_index = 0
        pool = 0

        # Get first amp number, fake msa_index
        self.assertEqual(1, multiplex.get_next_amplicon_number(100))
        self.assertEqual(1, multiplex.get_next_amplicon_number(0))

        # Create a primerpair
        primerpair = PrimerPair(
            FKmer([self.nCoV_2019_18_LEFT_0.encode()], 50),
            RKmer([b"TA"], 100),
            msa_index,
        )
        primerpair.pool = pool

        multiplex.add_primer_pair_to_pool(primerpair, pool, msa_index)

        # Check the amp num has updated
        self.assertEqual(2, multiplex.get_next_amplicon_number(0))

        # Add massive amp num
        primerpair2 = PrimerPair(FKmer([b"TA"], 50), RKmer([b"TA"], 100), 1)
        primerpair2.amplicon_number = 100
        multiplex.add_primer_pair_to_pool(primerpair2, pool, 1)

        self.assertEqual(101, multiplex.get_next_amplicon_number(1))


if __name__ == "__main__":
    unittest.main()
