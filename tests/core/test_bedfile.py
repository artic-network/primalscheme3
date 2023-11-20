import unittest
import pathlib
import hashlib
import os

from primalscheme3.core.bedfiles import (
    read_in_bedprimerpairs,
    BedPrimerPair,
    re_primer_name,
    read_in_bedlines,
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
    def test_round_trip(self):
        """
        Test that the bedfile can be read in and written out without any changes
        """

        input_path = pathlib.Path("tests/core/test_primer.bed").absolute()

        # Read in the bedfile
        bedprimerpairs: list[BedPrimerPair] = read_in_bedprimerpairs(input_path)
        # Hash the input bedfile
        input_hash = hashlib.file_digest(open(input_path, "rb"), "md5").hexdigest()

        # Create the bedstring
        bed_str = "".join([str(x) for x in bedprimerpairs]).strip()
        # Hash the output bedfile
        output_hash = hashlib.md5(bed_str.encode()).hexdigest()

        with open("tests/test_out.bed", "w") as f:
            f.write(bed_str)

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

        bedlines = read_in_bedlines(input_path)

        first_bedline = bedlines[0]
        # Check all primers were read in
        self.assertEqual(len(bedlines), 160)

        # Check the first line matches the expected output
        self.assertEqual(
            first_bedline.ref,
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


if __name__ == "__main__":
    unittest.main()
