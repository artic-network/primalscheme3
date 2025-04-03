import pathlib
import tempfile
import unittest

from primalscheme3.core.config import Config
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.panel.panel_main import panelcreate
from primalscheme3.scheme.scheme_main import schemecreate


class TestMain(unittest.TestCase):
    """
    This test the main functions of the primalscheme3 package.
    Doesn't validate outputs but checks if the functions run without errors.
    """

    msa_paths = [pathlib.Path("./tests/core/test_mismatch.fasta").absolute()]

    # Avoid using matchdb
    config = Config()
    config.use_matchdb = False

    def check_file(self, path):
        self.assertTrue(path.is_file())
        self.assertTrue(path.stat().st_size > 0)

    def test_schemecreate(self):
        with tempfile.TemporaryDirectory(
            dir="tests/integration", suffix="-schemecreate"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)

            # test_empty_path = tempdir_path / "empty"
            # test_empty_path.touch()
            # self.check_file(test_empty_path)

            # Run Scheme Create
            pm = ProgressManager()
            schemecreate(
                msa=self.msa_paths,
                output_dir=tempdir_path,
                pm=pm,
                config=self.config,
                force=True,
                offline_plots=False,
            )
            # Check for output files
            self.check_file(tempdir_path / "primer.bed")
            self.check_file(tempdir_path / "reference.fasta")
            self.check_file(tempdir_path / "plot.html")
            self.check_file(tempdir_path / "primer.html")
            self.check_file(tempdir_path / "config.json")

    def test_panelcreate_all(self):
        with tempfile.TemporaryDirectory(
            dir="tests/integration", suffix="-panelcreate"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)
            # Run Panel Create
            pm = ProgressManager()
            panelcreate(
                msa=self.msa_paths,
                output_dir=tempdir_path,
                config=self.config,
                pm=pm,
                force=True,
            )
            # Check for output files
            self.check_file(tempdir_path / "primer.bed")
            self.check_file(tempdir_path / "reference.fasta")
            self.check_file(tempdir_path / "plot.html")
            self.check_file(tempdir_path / "primer.html")
            self.check_file(tempdir_path / "config.json")
