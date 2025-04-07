import hashlib
import pathlib
import unittest

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
        input_hash = hashlib.file_digest(open(input_path, "rb"), "md5").hexdigest()

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
        input_hash = hashlib.file_digest(open(input_path, "rb"), "md5").hexdigest()

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
