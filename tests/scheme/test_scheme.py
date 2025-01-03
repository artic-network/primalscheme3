import pathlib
import unittest

from primalscheme3.core.classes import FKmer, MatchDB, PrimerPair, RKmer
from primalscheme3.core.config import Config
from primalscheme3.core.msa import MSA
from primalscheme3.scheme.classes import Scheme


class TestScheme(unittest.TestCase):
    db_path = pathlib.Path("./tests/core/mulitplex").absolute()
    matchdb = MatchDB(db_path, [], 30)  # Create an empty matchdb
    inputfile_path = pathlib.Path("./tests/core/test_mismatch.fasta").absolute()
    msa = MSA("test", inputfile_path, 0, "first", None)

    # Create a config dict
    config = Config()

    def test_get_leading_coverage_edge(self):
        """
        Test the method get_leading_coverage_edge
        """
        scheme = Scheme(
            config=self.config, matchDB=self.matchdb, msa_dict={0: self.msa}
        )
        primerpair = PrimerPair(FKmer(10, ["A"]), RKmer(20, ["T"]), 0)

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
        scheme = Scheme(
            config=self.config, matchDB=self.matchdb, msa_dict={0: self.msa}
        )
        primerpair = PrimerPair(FKmer(10, ["AA"]), RKmer(20, ["TT"]), 0)

        # Add a primerpair to pool 0
        scheme.add_primer_pair_to_pool(primerpair, 0, 0)

        # Check that the leading coverage edge is correct
        self.assertEqual(
            scheme.get_leading_amplicon_edge(),
            22,
        )

    def test_find_ol_primerpairs(self):
        """
        Test the method find_ol_primerpairs to produce the correct primerpairs
        """
        scheme = Scheme(
            config=self.config, matchDB=self.matchdb, msa_dict={0: self.msa}
        )
        primerpair = PrimerPair(FKmer(10, ["AA"]), RKmer(20, ["TT"]), 0)
        # Add a primerpair to pool 0
        scheme.add_primer_pair_to_pool(primerpair, 0, 0)

        # Create some overlapping primerpairs
        all_ol_primerpair = [
            PrimerPair(FKmer(x, ["AAA"]), RKmer(x + 100, ["TTT"]), 0)
            for x in range(50, 300, 10)
        ]
        # See which primers could ol
        pos_ol_primerpair = scheme.find_ol_primerpairs(
            all_ol_primerpair, self.config.min_overlap
        )

        # Make sure all primers have an overlap
        self.assertTrue(
            all(
                x.fprimer.end <= primerpair.rprimer.start - self.config.min_overlap
                for x in pos_ol_primerpair
            )
        )


if __name__ == "__main__":
    unittest.main()
