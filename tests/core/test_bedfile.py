import hashlib
import pathlib
import unittest

from primalscheme3.core.bedfiles import (
    re_primer_name,
    read_in_bedlines,
    read_in_bedprimerpairs,
)


class Test_RePrimerName(unittest.TestCase):
    def test_re_primer_name(self):
        # Test case 1: Valid input string
        string = "amplicon_1_RIGHT_1"
        expected_output = ["1", "RIGHT"]
        self.assertEqual(re_primer_name(string), expected_output)

        # Test case 2: Valid input string with leading and trailing spaces
        string = "  amplicon_2_LEFT_0  "
        expected_output = ["2", "LEFT"]
        self.assertEqual(re_primer_name(string), expected_output)

        # Test case 3: Invalid input string
        string = "invalid_string"
        expected_output = None
        self.assertEqual(re_primer_name(string), expected_output)

        # Test case 4: Empty input string
        string = ""
        expected_output = None
        self.assertEqual(re_primer_name(string), expected_output)


class Test_ReadInBedFile(unittest.TestCase):
    def test_round_trip_no_header(self):
        """
        Test that the bedfile can be read in and written out without any changes
        """

        input_path = pathlib.Path(
            "tests/test_data/test_bedfile_no_header.bed"
        ).absolute()

        # Read in the bedfile
        bedprimerpairs, _headers = read_in_bedprimerpairs(input_path)
        # Hash the input bedfile
        input_hash = hashlib.file_digest(open(input_path, "rb"), "md5").hexdigest()

        # Create the bedstring
        bed_list = []
        for headerline in _headers:
            if not headerline.startswith("#"):
                headerline = "# " + headerline
            bed_list.append(headerline.strip())

        for bedprimerpair in bedprimerpairs:
            bed_list.append(bedprimerpair.to_bed().strip())
        bed_str = "\n".join(bed_list)

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
        bedprimerpairs, _headers = read_in_bedprimerpairs(input_path)
        # Hash the input bedfile
        input_hash = hashlib.file_digest(open(input_path, "rb"), "md5").hexdigest()

        # Create the bedstring
        bed_list = []
        for headerline in _headers:
            if not headerline.startswith("#"):
                headerline = "# " + headerline
            bed_list.append(headerline.strip())

        for bedprimerpair in bedprimerpairs:
            bed_list.append(bedprimerpair.to_bed().strip())
        bed_str = "\n".join(bed_list)

        # Hash the output bedfile
        output_hash = hashlib.md5(bed_str.encode()).hexdigest()

        # with open("tests/test_data/test_bedfile_header_out.bed", "w") as outfile:
        #     outfile.write(bed_str)

        self.assertEqual(input_hash, output_hash)


class Test_ReadInBedlines(unittest.TestCase):
    def test_readin(self):
        """
        Test that a bedline can be read in correctly

        First bedline in test_primer.bed:
        MN908947.3	23	50	50dcc73b_0_LEFT_0	1	+	GTAACAAACCAACCAACTTTCGATCTC

        160 lines in test_primer.bed
        """
        input_path = pathlib.Path("tests/core/test_primer.bed").absolute()

        bedlines, _headers = read_in_bedlines(input_path)

        first_bedline = bedlines[0]
        # Check all primers were read in
        self.assertEqual(len(bedlines), 160)
        # Check empty headers
        self.assertEqual(_headers, [])

        # Check the first line matches the expected output
        self.assertEqual(
            first_bedline.chrom_name,
            "MN908947.3",
        )
        self.assertEqual(first_bedline.start, 23)
        self.assertEqual(first_bedline.end, 50)
        self.assertEqual(first_bedline.primername, "50dcc73b_0_LEFT_0")
        self.assertEqual(
            first_bedline.pool + 1, 1
        )  # Bedlines store the pool as zero indexed
        self.assertEqual(first_bedline.direction, "+")
        self.assertEqual(
            first_bedline.sequence,
            "GTAACAAACCAACCAACTTTCGATCTC",
        )

    def test_readin_header(self):
        """
        Test that a bedline with header can be read in correctly

        First bedline in test_primer.bed:
        #header
        # version 2.0 primer.bed
        MN908947.3	23	50	50dcc73b_0_LEFT_0	1	+	GTAACAAACCAACCAACTTTCGATCTC

        160 lines in test_primer.bed
        """
        input_path = pathlib.Path("tests/core/test_primer_header.bed").absolute()

        bedlines, headers = read_in_bedlines(input_path)

        # Check all primers were read in
        self.assertEqual(headers, ["#header", "# version 2.0 primer.bed"])

        # Check the first line matches the expected output
        self.assertEqual(len(bedlines), 160)


if __name__ == "__main__":
    unittest.main()
