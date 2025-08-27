import hashlib
import pathlib
import unittest

from primalschemers import FKmer, RKmer  # type: ignore

from primalscheme3.core.bedfiles import (
    read_bedlines_to_bedprimerpairs,
)


class Test_ReadInBedFile(unittest.TestCase):
    def test_round_trip_no_header(self):
        """
        Test that the bedfile can be read in and written out without any changes
        """

        input_path = pathlib.Path(
            "tests/test_data/test_bedfile_no_header.bed"
        ).absolute()

        # Read in the bedfile
        bedprimerpairs, _headers = read_bedlines_to_bedprimerpairs(input_path)
        # Hash the input bedfile
        with open(input_path, "rb") as file:
            input_hash = hashlib.file_digest(file, "md5").hexdigest()

        # Create the bedstring
        bed_list = []
        for headerline in _headers:
            if not headerline.startswith("#"):
                headerline = "# " + headerline
            bed_list.append(headerline.strip())

        for bedprimerpair in bedprimerpairs:
            bed_list.append(bedprimerpair.to_bed().strip())
        bed_str = "\n".join(bed_list) + "\n"

        # Hash the output bedfile
        output_hash = hashlib.md5(bed_str.encode()).hexdigest()

        # with open("tests/test_data/test_bedfile_No_header_out.bed", "w") as outfile:
        #   outfile.write(bed_str)

        self.assertEqual(input_hash, output_hash)

    def test_round_trip_header(self):
        """
        Test that the bedfile can be read in and written out without any changes
        """

        input_path = pathlib.Path("tests/test_data/test_bedfile_header.bed").absolute()

        # Read in the bedfile
        bedprimerpairs, _headers = read_bedlines_to_bedprimerpairs(input_path)
        # Hash the input bedfile
        with open(input_path, "rb") as file:
            input_hash = hashlib.file_digest(file, "md5").hexdigest()

        # Create the bedstring
        bed_list = []
        for headerline in _headers:
            if not headerline.startswith("#"):
                headerline = "# " + headerline
            bed_list.append(headerline.strip())

        for bedprimerpair in bedprimerpairs:
            bed_list.append(bedprimerpair.to_bed().strip())
        bed_str = "\n".join(bed_list) + "\n"

        # Hash the output bedfile
        output_hash = hashlib.md5(bed_str.encode()).hexdigest()

        # with open("tests/test_data/test_bedfile_header_out.bed", "w") as outfile:
        #     outfile.write(bed_str)

        self.assertEqual(input_hash, output_hash)

    def test_kmer_counts(self):
        """Test that the count attr is only set of KMers with counts"""
        chrom = "chrom"
        amplicon_prefix = "amp_5"
        pool = 1

        # With count
        fkmer = FKmer([b"CGATCGAC"], 20, counts=[10])
        self.assertEqual(
            fkmer.to_bed(chrom=chrom, amplicon_prefix=amplicon_prefix, pool=pool),
            f"{chrom}\t12\t20\t{amplicon_prefix}_LEFT_1\t{pool}\t+\tCGATCGAC\tpc=10\n",
        )
        # Without count
        fkmer = FKmer([b"CGATCGAC"], 20)
        self.assertEqual(
            fkmer.to_bed(chrom=chrom, amplicon_prefix=amplicon_prefix, pool=pool),
            f"{chrom}\t12\t20\t{amplicon_prefix}_LEFT_1\t{pool}\t+\tCGATCGAC\n",
        )

        # With count
        rkmer = RKmer([b"CGATCGAC"], 20, counts=[10])
        self.assertEqual(
            rkmer.to_bed(chrom=chrom, amplicon_prefix=amplicon_prefix, pool=pool),
            f"{chrom}\t20\t28\t{amplicon_prefix}_RIGHT_1\t{pool}\t-\tCGATCGAC\tpc=10\n",
        )
        # Without count
        rkmer = RKmer([b"CGATCGAC"], 20)
        self.assertEqual(
            rkmer.to_bed(chrom=chrom, amplicon_prefix=amplicon_prefix, pool=pool),
            f"{chrom}\t20\t28\t{amplicon_prefix}_RIGHT_1\t{pool}\t-\tCGATCGAC\n",
        )


if __name__ == "__main__":
    unittest.main()
