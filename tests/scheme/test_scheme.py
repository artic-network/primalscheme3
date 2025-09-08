import pathlib
import unittest

from primalschemers import FKmer, RKmer  # type: ignore

from primalscheme3.core.classes import PrimerPair
from primalscheme3.core.config import Config
from primalscheme3.core.mismatches import MatchDB
from primalscheme3.core.msa import MSA
from primalscheme3.scheme.classes import Scheme, SchemeReturn


class TestScheme(unittest.TestCase):
    def setUp(self) -> None:
        self.config = Config()
        self.db_path = pathlib.Path("./tests/core/mulitplex").absolute()
        self.matchdb = MatchDB(self.db_path, [], self.config)  # Create an empty matchdb
        self.inputfile_path = pathlib.Path(
            "./tests/core/test_mismatch.fasta"
        ).absolute()

        self.msa = MSA("test", self.inputfile_path, 0, None, self.config)

        return super().setUp()

    def test_get_leading_coverage_edge(self):
        """
        Test the method get_leading_coverage_edge
        """
        scheme = Scheme(
            config=self.config, matchDB=self.matchdb, msa_dict={0: self.msa}
        )
        primerpair = PrimerPair(
            FKmer([b"A"], 10),
            RKmer([b"T"], 20),
            0,
        )

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
        primerpair = PrimerPair(FKmer([b"AA"], 10), RKmer([b"TT"], 20), 0)

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
        primerpair = PrimerPair(FKmer([b"AA"], 10), RKmer([b"TT"], 20), 0)
        # Add a primerpair to pool 0
        scheme.add_primer_pair_to_pool(primerpair, 0, 0)

        # Create some overlapping primerpairs
        all_ol_primerpair = [
            PrimerPair(FKmer([b"AAA"], x), RKmer([b"TTT"], x + 100), 0)
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

    def test_backtrack(self):
        self.config.min_overlap = 0
        scheme = Scheme(
            config=self.config, matchDB=self.matchdb, msa_dict={0: self.msa}
        )
        primerpair = PrimerPair(
            FKmer([b"ACCAACGATGGTGTGTCCAT"], 10),
            RKmer([b"CTTGTCGAACCGCATACCCT"], 50),
            0,
        )
        all_ol_primerpair = [
            PrimerPair(
                FKmer([b"GCGACGGGTACGAGTGGTCT"], x),
                RKmer([b"CGTTCCATTGCATCGCGATCTC"], x + 100),
                0,
            )
            for x in range(40, 300, 10)
        ]
        # Add the first primerpair
        scheme.add_primer_pair_to_pool(primerpair, 0, 0)

        # Add second blocking primerpair.
        # Same right end prevents ol
        block_pp = PrimerPair(FKmer([b"AA"], 15), RKmer([b"TT"], 50), 0)
        scheme.add_primer_pair_to_pool(block_pp, 1, 0)

        # Show no ol primerpairs can be added
        self.assertEqual(
            scheme.try_ol_primerpairs(all_ol_primerpair, 0),
            SchemeReturn.NO_OL_PRIMERPAIR,
        )

        # Show bt can solve it
        self.assertEqual(
            scheme.try_backtrack(all_ol_primerpair, 0), SchemeReturn.ADDED_BACKTRACKED
        )
        # Blocking as been replaced with new pp
        self.assertTrue(scheme._last_pp_added[-1] != block_pp)
        self.assertTrue(scheme._last_pp_added[-1] in all_ol_primerpair)


if __name__ == "__main__":
    unittest.main()
