import unittest

from primalscheme3.core.classes import FKmer, RKmer, PrimerPair


class Test_FKmer(unittest.TestCase):
    def test_creation(self):
        seqs = {"ATGC", "ATGC", "ATGC", "ACACAA"}
        end = 100

        # Test case 1: Valid input
        fkmer = FKmer(end=end, seqs=seqs)

        # Test asignments
        self.assertEqual(fkmer.seqs, seqs)
        self.assertEqual(fkmer.end, end)
        self.assertEqual(fkmer._starts, {end - 4, end - 6})

    def test_len(self):
        seqs = {"ATGC", "ATGC", "ATGC"}
        end = 100

        # Test case 1: Valid input
        fkmer = FKmer(end=end, seqs=seqs)

        # Test asignments
        self.assertEqual(fkmer.len(), {4})

    def test_string_single(self):
        seqs = {"ATGC"}
        end = 100
        reference = "reference"
        amplicon_prefix = "amplicon_prefix"
        pool = "pool"

        # Test case 1: Valid input
        fkmer = FKmer(end=end, seqs=seqs)

        # Test asignments
        self.assertEqual(
            fkmer.__str__(reference, amplicon_prefix, pool),
            "reference\t96\t100\tamplicon_prefix_LEFT_0\tpool\t+\tATGC\n",
        )

    def test_string_multiple(self):
        seqs = {"ATGC", "ATGCA"}
        end = 100
        reference = "reference"
        amplicon_prefix = "amplicon_prefix"
        pool = "pool"

        # Test case 1: Valid input
        fkmer = FKmer(end=end, seqs=seqs)

        # Test asignments
        self.assertEqual(
            fkmer.__str__(reference, amplicon_prefix, pool),
            "reference\t96\t100\tamplicon_prefix_LEFT_0\tpool\t+\tATGC\nreference\t95\t100\tamplicon_prefix_LEFT_1\tpool\t+\tATGCA\n",
        )


class Test_RKmer(unittest.TestCase):
    def test_create(self):
        seqs = {"ATGC", "ATGC", "ATGC"}
        start = 100

        # Test case 1: Valid input
        rkmer = RKmer(start=start, seqs=seqs)

        # Test asignments
        self.assertEqual(rkmer.seqs, seqs)
        self.assertEqual(rkmer.start, start)
        self.assertEqual(rkmer._ends, {start + 4})

    def test_len(self):
        seqs = {"ATGC", "ATGC", "ATGC"}
        start = 100

        # Test case 1: Valid input
        rkmer = RKmer(start=start, seqs=seqs)

        # Test asignments
        self.assertEqual(rkmer.len(), {4})

    def test_ends(self):
        seqs = {"ATGC", "ATGC", "ATGCAA"}
        start = 100

        # Test case 1: Valid input
        rkmer = RKmer(start=start, seqs=seqs)

        # Test asignments
        self.assertEqual(rkmer.ends(), {104, 106})

    def test_string_single(self):
        seqs = {"ATGC"}
        start = 100
        reference = "reference"
        amplicon_prefix = "amplicon_prefix"
        pool = "pool"

        # Test case 1: Valid input
        rkmer = RKmer(start=start, seqs=seqs)

        # Test asignments
        self.assertEqual(
            rkmer.__str__(reference, amplicon_prefix, pool),
            "reference\t100\t104\tamplicon_prefix_RIGHT_0\tpool\t-\tATGC\n",
        )

    def test_string_multiple(self):
        seqs = {"ATGC", "ATGCA"}
        start = 100
        reference = "reference"
        amplicon_prefix = "amplicon_prefix"
        pool = "pool"

        # Test case 1: Valid input
        rkmer = RKmer(start=start, seqs=seqs)

        # Test asignments
        self.assertEqual(
            rkmer.__str__(reference, amplicon_prefix, pool),
            "reference\t100\t104\tamplicon_prefix_RIGHT_0\tpool\t-\tATGC\nreference\t100\t105\tamplicon_prefix_RIGHT_1\tpool\t-\tATGCA\n",
        )


class Test_PrimerPair(unittest.TestCase):
    def test_create(self):
        fkmer = FKmer(end=100, seqs={"ATGC"})
        rkmer = RKmer(start=1000, seqs={"ATGC"})
        msa_index = 0

        # Test case 1: Valid input
        primerpair = PrimerPair(fprimer=fkmer, rprimer=rkmer, msa_index=msa_index)

        # Test asignments
        self.assertEqual(primerpair.fprimer, fkmer)
        self.assertEqual(primerpair.rprimer, rkmer)
        self.assertEqual(primerpair.msa_index, msa_index)

    def test_set_amplicon_number(self):
        fkmer = FKmer(end=100, seqs={"ATGC"})
        rkmer = RKmer(start=1000, seqs={"ATGC"})
        msa_index = 0

        # Test case 1: Valid input
        primerpair = PrimerPair(fprimer=fkmer, rprimer=rkmer, msa_index=msa_index)
        primerpair.set_amplicon_number(1)

        # Test asignments
        self.assertEqual(primerpair.amplicon_number, 1)

    def test_all_seqs(self):
        fkmer = FKmer(end=100, seqs={"ACTAGCTAGCTAGCA"})
        rkmer = RKmer(start=1000, seqs={"ATCGATCGGTAC"})
        msa_index = 0

        # Test case 1: Valid input
        primerpair = PrimerPair(fprimer=fkmer, rprimer=rkmer, msa_index=msa_index)

        # Test asignments
        self.assertEqual(primerpair.all_seqs(), ["ACTAGCTAGCTAGCA", "ATCGATCGGTAC"])

    def test_to_bed(self):
        fkmer = FKmer(end=100, seqs={"ACTAGCTAGCTAGCA"})
        rkmer = RKmer(start=1000, seqs={"ATCGATCGGTAC"})
        msa_index = 0

        # Test case 1: Valid input
        primerpair = PrimerPair(fprimer=fkmer, rprimer=rkmer, msa_index=msa_index)
        primerpair.pool = 0
        primerpair.set_amplicon_number(0)

        # Test asignments
        expected_pool = primerpair.pool + 1
        expected_refname = "reference"
        expected_amplicon_prefix = "amplicon"

        primerpair.chrom_name = expected_refname
        primerpair.amplicon_prefix = expected_amplicon_prefix

        expected_str = f"{expected_refname}\t85\t100\t{expected_amplicon_prefix}_0_LEFT_0\t{expected_pool}\t+\tACTAGCTAGCTAGCA\n{expected_refname}\t1000\t1012\t{expected_amplicon_prefix}_0_RIGHT_0\t{expected_pool}\t-\tATCGATCGGTAC\n"

        self.assertEqual(primerpair.to_bed(), expected_str)


if __name__ == "__main__":
    unittest.main()
