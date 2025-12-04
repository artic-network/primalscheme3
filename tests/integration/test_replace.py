import json
import pathlib
import tempfile
import unittest

from primalbedtools.scheme import Scheme

from primalscheme3.core.config import Config
from primalscheme3.replace.replace import ReplaceRunModes, replace


class TestRepairMode(unittest.TestCase):
    input_bedfile = pathlib.Path("tests/test_data/test_scheme/primer.bed").absolute()
    input_config = pathlib.Path("tests/test_data/test_scheme/config.json").absolute()
    input_msa = pathlib.Path(
        "tests/test_data/test_scheme/work/reference.fasta"
    ).absolute()

    def setUp(self) -> None:
        # Read in the config file
        with open(self.input_config) as file:
            _cfg: dict = json.load(file)
        self.config = Config(**_cfg)
        return super().setUp()

    def test_ol_replace(self):
        with tempfile.TemporaryDirectory(
            dir="tests/integration", suffix="-replace"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)

            self.config.amplicon_size_max = 900

            # Run the program
            replace(
                config=self.config,
                primerbed=self.input_bedfile,
                primername="1540b43d_11_RIGHT_1",
                msapath=self.input_msa,
                pm=None,
                output=tempdir_path,
                force=True,
                mode=ReplaceRunModes.AddBest,
            )

            # Read in the output primerbed file
            out_primer_bed = tempdir_path / "primer.bed"
            self.assertTrue(out_primer_bed.exists())

            out_scheme = Scheme.from_file(str(tempdir_path / "primer.bed"))
            for bl in out_scheme.bedlines:
                self.assertFalse(bl.amplicon_name == "1540b43d_11")


if __name__ == "__main__":
    unittest.main()
