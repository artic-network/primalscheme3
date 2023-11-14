import unittest
import pathlib
import hashlib

from primalscheme3.core.bedfiles import read_in_bedfile, BedPrimerPair


class Test_ReadInBedFile(unittest.TestCase):
    def test_round_trip(self):
        """
        Test that the bedfile can be read in and written out without any changes
        """
        import os

        print(os.getcwd())

        input_path = pathlib.Path("tests/test_primer.bed").absolute()

        # Read in the bedfile
        bedprimerpairs: list[BedPrimerPair] = read_in_bedfile(input_path)
        # Hash the input bedfile
        input_hash = hashlib.file_digest(open(input_path, "rb"), "md5").hexdigest()

        # Create the bedstring
        bed_str = "".join([str(x) for x in bedprimerpairs]).strip()
        # Hash the output bedfile
        output_hash = hashlib.md5(bed_str.encode()).hexdigest()

        with open("tests/test_out.bed", "w") as f:
            f.write(bed_str)

        print(bed_str)

        self.assertEqual(input_hash, output_hash)


if __name__ == "__main__":
    unittest.main()
