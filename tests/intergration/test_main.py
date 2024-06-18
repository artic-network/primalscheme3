import pathlib
import unittest

from primalscheme3.core.config import config_dict as cfg
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.scheme.scheme_main import schemecreate


class TestMain(unittest.TestCase):
    msa_paths = [pathlib.Path("./tests/core/test_mismatch.fasta").absolute()]
    outpath = pathlib.Path("./tests/intergration/output").absolute()
    pm = ProgressManager()

    def setUp(self) -> None:
        self.outpath.mkdir(exist_ok=True, parents=True)
        return super().setUp()

    def test_schemecreate(self):
        # Run Scheme Create
        schemecreate(
            argmsa=self.msa_paths,
            ampliconsizemax=440,
            ampliconsizemin=360,
            dimerscore=cfg["dimerscore"],
            output_dir=self.outpath,
            primer_gc_max=cfg["primer_gc_max"],
            primer_gc_min=cfg["primer_gc_min"],
            primer_tm_max=cfg["primer_tm_max"],
            primer_tm_min=cfg["primer_tm_min"],
            ncores=1,
            minoverlap=cfg["minoverlap"],
            npools=2,
            reducekmers=False,
            minbasefreq=0,
            circular=False,
            backtrack=False,
            ignore_n=False,
            pm=self.pm,
            force=True,
        )
