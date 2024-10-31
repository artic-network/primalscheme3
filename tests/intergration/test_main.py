import pathlib
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
    outpath = pathlib.Path("./tests/integration/output").absolute()

    # Avoid using matchdb
    config = Config()
    config.use_matchdb = False

    def setUp(self) -> None:
        self.outpath.mkdir(exist_ok=True, parents=True)
        return super().setUp()

    def test_schemecreate(self):
        # Run Scheme Create
        pm = ProgressManager()
        schemecreate(
            msa=self.msa_paths,
            output_dir=self.outpath / "schemecreate",
            pm=pm,
            config=self.config,
            force=True,
            offline_plots=False,
        )

    def test_panelcreate_all(self):
        # Run Panel Create
        pm = ProgressManager()
        panelcreate(
            msa=self.msa_paths,
            output_dir=self.outpath / "panelcreate_all",
            config=self.config,
            pm=pm,
            force=True,
        )
