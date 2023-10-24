import unittest
import pathlib

from primal_digest.classes import Scheme, MatchDB, PrimerPair, FKmer, RKmer
import primal_digest.config as config


class TestScheme(unittest.TestCase):
    db_path = pathlib.Path("tests").absolute()
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
        scheme = Scheme(cfg=self.cfg, matchDB=self.matchdb)
        current_pool = scheme._current_pool
        next_pool = scheme.next_pool()
        self.assertEqual(
            current_pool + 1,
            next_pool,
        )

    def test_next_pool_1_pool(self):
        """
        Test if calling the next_pool method returns the correct pool
        """
        self.cfg["npools"] = 1
        scheme = Scheme(cfg=self.cfg, matchDB=self.matchdb)
        current_pool = scheme._current_pool
        next_pool = scheme.next_pool()
        self.assertEqual(
            current_pool,
            next_pool,
        )

    def test_add_primer_pair_to_pool(self):
        """
        Test if method add_primer_pair_to_pool does whats expected
        """
        self.cfg["npools"] = 2
        scheme = Scheme(cfg=self.cfg, matchDB=self.matchdb)
        pp_msa_index = 0

        primerpair = PrimerPair(FKmer(10, "A"), RKmer(20, "T"), None)
        # Add a primerpair to pool 0
        scheme.add_primer_pair_to_pool(primerpair, 0, pp_msa_index)

        # Check that the primerpair has beed added to _last_pp_added
        self.assertEqual(
            scheme._last_pp_added[-1],
            primerpair,
        )
        # Check that the primerpair has been added to the correct pool
        self.assertEqual(scheme._pools[0], [primerpair])
        # Check that current pool has updated
        self.assertEqual(scheme._current_pool, 1)
        # Check that the primerpair has had its msa_index updated
        self.assertEqual(primerpair.msa_index, pp_msa_index)
        # Check amplicon number has been assigned
        self.assertEqual(scheme._last_pp_added[-1].amplicon_number, 0)

    def test_get_leading_coverage_edge(self):
        """
        Test the method get_leading_coverage_edge
        """
        self.cfg["npools"] = 2
        scheme = Scheme(cfg=self.cfg, matchDB=self.matchdb)
        primerpair = PrimerPair(FKmer(10, "A"), RKmer(20, "T"), None)

        # Add a primerpair to pool 0
        scheme.add_primer_pair_to_pool(primerpair, 0, 0)

        # Check that the leading coverage edge is correct
        self.assertEqual(
            scheme.get_leading_coverage_edge(),
            20,
        )

    def test_get_leading_amplicon_edge(self):
        """
        Test the method get_leading_coverage_edge
        """
        self.cfg["npools"] = 2
        scheme = Scheme(cfg=self.cfg, matchDB=self.matchdb)
        primerpair = PrimerPair(FKmer(10, "AA"), RKmer(20, "TT"), None)

        # Add a primerpair to pool 0
        scheme.add_primer_pair_to_pool(primerpair, 0, 0)

        # Check that the leading coverage edge is correct
        self.assertEqual(
            scheme.get_leading_amplicon_edge(),
            21,
        )

    def test_find_ol_primerpairs(self):
        """
        Test the method find_ol_primerpairs to produce the correct primerpairs
        """
        self.cfg["npools"] = 2
        scheme = Scheme(cfg=self.cfg, matchDB=self.matchdb)
        primerpair = PrimerPair(FKmer(100, "AA"), RKmer(200, "TT"), None)

        # Add a primerpair to pool 0
        scheme.add_primer_pair_to_pool(primerpair, 0, 0)

        # Create some overlapping primerpairs
        all_ol_primerpair = [
            PrimerPair(FKmer(x, "AAA"), RKmer(x + 100, "TTT"), None)
            for x in range(50, 300, 10)
        ]
        # See which primers could ol
        pos_ol_primerpair = scheme.find_ol_primerpairs(all_ol_primerpair)

        # Make sure all primers have an overlap
        self.assertTrue(
            all(
                (
                    x.fprimer.end <= primerpair.rprimer.start - self.cfg["min_overlap"]
                    for x in pos_ol_primerpair
                )
            )
        )


if __name__ == "__main__":
    unittest.main()
