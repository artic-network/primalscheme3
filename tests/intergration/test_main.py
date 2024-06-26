import pathlib
import unittest

from primalscheme3.core.config import Config
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.scheme.scheme_main import schemecreate


class TestMain(unittest.TestCase):
    msa_paths = [pathlib.Path("./tests/core/test_mismatch.fasta").absolute()]
    outpath = pathlib.Path("./tests/intergration/output").absolute()
    pm = ProgressManager()
    config = Config()

    def setUp(self) -> None:
        self.outpath.mkdir(exist_ok=True, parents=True)
        return super().setUp()

    def test_schemecreate(self):
        # Run Scheme Create
        schemecreate(
            msa=self.msa_paths,
            output_dir=self.outpath,
            pm=self.pm,
            config=self.config,
            force=True,
        )
