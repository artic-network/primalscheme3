import unittest

from primalschemers import FKmer, RKmer  # type: ignore

from primalscheme3.core.classes import PrimerPair


class Test_PrimerPair(unittest.TestCase):
    def test_create(self):
        fkmer = FKmer([b"ATGC"], 100)
        rkmer = RKmer([b"ATGC"], 1000)
        msa_index = 0

        # Test case 1: Valid input
        primerpair = PrimerPair(fprimer=fkmer, rprimer=rkmer, msa_index=msa_index)

        # Test asignments
        self.assertEqual(primerpair.fprimer, fkmer)
        self.assertEqual(primerpair.rprimer, rkmer)
        self.assertEqual(primerpair.msa_index, msa_index)

    def test_set_amplicon_number(self):
        fkmer = FKmer([b"ATGC"], 100)
        rkmer = RKmer([b"ATGC"], 1000)
        msa_index = 0

        # Test case 1: Valid input
        primerpair = PrimerPair(fprimer=fkmer, rprimer=rkmer, msa_index=msa_index)
        primerpair.set_amplicon_number(1)

        # Test asignments
        self.assertEqual(primerpair.amplicon_number, 1)

    def test_all_seqs(self):
        fkmer = FKmer([b"ACTAGCTAGCTAGCA"], 100)
        rkmer = RKmer([b"ATCGATCGGTAC"], 1000)
        msa_index = 0

        primerpair = PrimerPair(fprimer=fkmer, rprimer=rkmer, msa_index=msa_index)

        # Test asignments
        self.assertEqual(primerpair.all_seqs(), ["ACTAGCTAGCTAGCA", "ATCGATCGGTAC"])

    def test_all_seqs_bytes(self):
        fkmer = FKmer([b"ACTAGCTAGCTAGCA"], 100)
        rkmer = RKmer([b"ATCGATCGGTAC"], 1000)
        msa_index = 0

        primerpair = PrimerPair(fprimer=fkmer, rprimer=rkmer, msa_index=msa_index)

        # Test asignments
        seqs = primerpair.all_seqs()
        seqs_bytes = primerpair.all_seq_bytes()

        for s, sb in zip(seqs, seqs_bytes, strict=False):
            self.assertEqual(s.encode(), sb)

    def test_to_bed(self):
        fkmer = FKmer([b"ACTAGCTAGCTAGCA"], 100)
        rkmer = RKmer([b"ATCGATCGGTAC"], 1000)
        msa_index = 0

        # Test case 1: Valid input
        primerpair = PrimerPair(fprimer=fkmer, rprimer=rkmer, msa_index=msa_index)
        primerpair.pool = 0
        primerpair.set_amplicon_number(0)

        # Test alignments
        expected_pool = primerpair.pool + 1
        expected_refname = "reference"
        expected_amplicon_prefix = "amplicon"

        primerpair.chrom_name = expected_refname
        primerpair.amplicon_prefix = expected_amplicon_prefix

        expected_str = f"{expected_refname}\t85\t100\t{expected_amplicon_prefix}_0_LEFT_1\t{expected_pool}\t+\tACTAGCTAGCTAGCA\n{expected_refname}\t1000\t1012\t{expected_amplicon_prefix}_0_RIGHT_1\t{expected_pool}\t-\tATCGATCGGTAC\n"

        self.assertEqual(primerpair.to_bed(), expected_str)


if __name__ == "__main__":
    unittest.main()
