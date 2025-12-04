import pathlib
import tempfile
import unittest

from primalscheme3.core.config import Config, MappingType
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.panel.panel_main import PanelRunModes, panelcreate
from primalscheme3.scheme.scheme_main import schemecreate


class TestMain(unittest.TestCase):
    """
    This test the main functions of the primalscheme3 package.
    Doesn't validate outputs but checks if the functions run without errors.
    """

    msa_paths = [pathlib.Path("./tests/core/test_mismatch.fasta").absolute()]
    region_path = pathlib.Path("./tests/core/test_mismatch.fasta.bed").absolute()
    input_bedfile = pathlib.Path("tests/test_data/test_scheme/primer.bed").absolute()
    input_config = pathlib.Path("tests/test_data/test_scheme/config.json").absolute()
    input_msa_path = pathlib.Path(
        "tests/test_data/test_scheme/work/reference.fasta"
    ).absolute()

    def setUp(self) -> None:
        self.config = Config()
        self.config.use_matchdb = False
        return super().setUp()

    def check_file(self, path):
        self.assertTrue(path.is_file())
        self.assertTrue(path.stat().st_size > 0)

    def test_schemecreate_first(self):
        with tempfile.TemporaryDirectory(
            dir="tests/integration", suffix="-schemecreate"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)

            # Modify config
            self.config.mapping = MappingType.FIRST

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

    def test_schemecreate_circular(self):
        with tempfile.TemporaryDirectory(
            dir="tests/integration", suffix="-schemecreate"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)

            # Modify config
            self.config.mapping = MappingType.FIRST
            self.config.circular = True

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

    def test_schemecreate_first_input_bed(self):
        with tempfile.TemporaryDirectory(
            dir="tests/integration", suffix="-schemecreate"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)

            # Modify config
            self.config.mapping = MappingType.FIRST

            # Run Scheme Create
            pm = ProgressManager()
            schemecreate(
                msa=self.msa_paths,
                output_dir=tempdir_path,
                pm=pm,
                config=self.config,
                force=True,
                offline_plots=False,
                input_bedfile=self.input_bedfile,
            )
            # Check for output files
            self.check_file(tempdir_path / "primer.bed")
            self.check_file(tempdir_path / "reference.fasta")
            self.check_file(tempdir_path / "plot.html")
            self.check_file(tempdir_path / "primer.html")
            self.check_file(tempdir_path / "config.json")

    def test_schemecreate_first_backtrack(self):
        with tempfile.TemporaryDirectory(
            dir="tests/integration", suffix="-schemecreate"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)

            # Modify config
            self.config.mapping = MappingType.FIRST
            self.config.backtrack = True

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

    def test_schemecreate_consensus(self):
        with tempfile.TemporaryDirectory(
            dir="tests/integration", suffix="-schemecreate"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)

            # Modify config
            self.config.mapping = MappingType.CONSENSUS

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

    def test_panelcreate_consensus_entropy(self):
        with tempfile.TemporaryDirectory(
            dir="tests/integration", suffix="-panelcreate"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)

            # Modify config
            self.config.mapping = MappingType.CONSENSUS
            mode = PanelRunModes.ENTROPY
            # Run Panel Create
            pm = ProgressManager()
            panelcreate(
                msa=self.msa_paths,
                output_dir=tempdir_path,
                config=self.config,
                pm=pm,
                force=True,
                mode=mode,
            )
            # Check for output files
            self.check_file(tempdir_path / "primer.bed")
            self.check_file(tempdir_path / "reference.fasta")
            self.check_file(tempdir_path / "plot.html")
            self.check_file(tempdir_path / "primer.html")
            self.check_file(tempdir_path / "config.json")

    def test_panelcreate_consensus_equal(self):
        with tempfile.TemporaryDirectory(
            dir="tests/integration", suffix="-panelcreate"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)

            # Modify config
            self.config.mapping = MappingType.CONSENSUS
            mode = PanelRunModes.EQUAL
            # Run Panel Create
            pm = ProgressManager()
            panelcreate(
                msa=self.msa_paths,
                output_dir=tempdir_path,
                config=self.config,
                pm=pm,
                force=True,
                mode=mode,
            )
            # Check for output files
            self.check_file(tempdir_path / "primer.bed")
            self.check_file(tempdir_path / "reference.fasta")
            self.check_file(tempdir_path / "plot.html")
            self.check_file(tempdir_path / "primer.html")
            self.check_file(tempdir_path / "config.json")

    def test_panelcreate_first_region_only(self):
        with tempfile.TemporaryDirectory(
            dir="tests/integration", suffix="-panelcreate"
        ) as tempdir:
            tempdir_path = pathlib.Path(tempdir)

            # Modify config
            self.config.mapping = MappingType.FIRST
            mode = PanelRunModes.REGION_ONLY
            # Run Panel Create
            pm = ProgressManager()
            panelcreate(
                msa=self.msa_paths,
                output_dir=tempdir_path,
                config=self.config,
                pm=pm,
                force=True,
                mode=mode,
                region_bedfile=self.region_path,
            )
            # Check for output files
            self.check_file(tempdir_path / "primer.bed")
            self.check_file(tempdir_path / "reference.fasta")
            self.check_file(tempdir_path / "plot.html")
            self.check_file(tempdir_path / "primer.html")
            self.check_file(tempdir_path / "config.json")
