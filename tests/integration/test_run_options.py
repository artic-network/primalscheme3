import json
import pathlib
import tempfile
import unittest

from primalscheme3.core.config import Config, MappingType
from primalscheme3.core.errors import UsageError
from primalscheme3.core.msa import parse_msa
from primalscheme3.core.primer_visual import primer_dct_heatmap
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.panel.panel_main import PanelRunModes, panelcreate
from primalscheme3.replace.replace import ReplaceRunModes, replace
from primalscheme3.scheme.scheme_main import schemecreate

MSA = pathlib.Path("tests/core/test_mismatch.fasta").absolute()
SCHEME_DIR = pathlib.Path("tests/test_data/test_scheme").absolute()


def amplicons(out: pathlib.Path) -> list[list[str]]:
    """The amplicon.bed rows: chrom, start, end, name, pool"""
    return [
        line.split("\t")
        for line in (out / "amplicon.bed").read_text().splitlines()
        if line and not line.startswith("#")
    ]


def renamed_copy(src: pathlib.Path, dest: pathlib.Path, prefix: str) -> pathlib.Path:
    """The MSA with each sequence id prefixed, so it reads as another MSA"""
    dest.write_text(src.read_text().replace(">", f">{prefix}"))
    return dest


def run_panel(out: pathlib.Path, config: Config, **kwargs) -> dict:
    panelcreate(
        output_dir=out,
        config=config,
        pm=ProgressManager(),
        force=True,
        offline_plots=False,
        **{"msa": [MSA], "mode": PanelRunModes.ENTROPY, **kwargs},
    )
    return json.loads((out / "config.json").read_text())


class TestSchemeAnnealing(unittest.TestCase):
    def test_use_annealing(self):
        config = Config(use_matchdb=False, use_annealing=True)
        with tempfile.TemporaryDirectory(dir="tests/integration") as d:
            out = pathlib.Path(d)
            schemecreate(
                msa=[MSA],
                output_dir=out,
                pm=ProgressManager(),
                config=config,
                force=True,
                offline_plots=False,
            )
            written = json.loads((out / "config.json").read_text())
            self.assertIsNotNone(written["primer_annealing_prop"])
            self.assertTrue(amplicons(out))


class TestPanelOptions(unittest.TestCase):
    def config(self, **kwargs) -> Config:
        return Config(amplicon_size=200, **kwargs)

    def test_mispriming_check(self):
        with tempfile.TemporaryDirectory(dir="tests/integration") as d:
            out = pathlib.Path(d)
            written = run_panel(out, self.config(use_matchdb=True))
            self.assertTrue(amplicons(out))
            self.assertGreater(written["pair_checks"]["OK"], 0)

    def test_max_amplicons(self):
        with tempfile.TemporaryDirectory(dir="tests/integration") as d:
            out = pathlib.Path(d)
            run_panel(out, self.config(use_matchdb=False), max_amplicons=2)
            self.assertEqual(len(amplicons(out)), 2)

    def test_max_amplicons_msa(self):
        with tempfile.TemporaryDirectory(dir="tests/integration") as d:
            d = pathlib.Path(d)
            copy = renamed_copy(MSA, d / "copy.fasta", "copy_")
            out = d / "out"
            run_panel(
                out,
                self.config(use_matchdb=False),
                msa=[MSA, copy],
                max_amplicons_msa=1,
            )
            per_chrom = {}
            for chrom, *_ in amplicons(out):
                per_chrom[chrom] = per_chrom.get(chrom, 0) + 1
            self.assertEqual(per_chrom, {"sequence1": 1, "copy_sequence1": 1})

    def test_max_amplicons_region_group(self):
        groups = {"a": [(100, 110), (300, 310)], "b": [(600, 610), (800, 810)]}
        with tempfile.TemporaryDirectory(dir="tests/integration") as d:
            d = pathlib.Path(d)
            regions = d / "regions.bed"
            regions.write_text(
                "".join(
                    f"sequence1\t{start}\t{end}\t{group}{i}\t10\t{group}\n"
                    for group, sites in groups.items()
                    for i, (start, end) in enumerate(sites)
                )
            )
            for limit, expected_max in ((None, 4), (1, 2)):
                with self.subTest(limit=limit):
                    out = d / f"out_{limit}"
                    run_panel(
                        out,
                        self.config(use_matchdb=False),
                        mode=PanelRunModes.REGION_ONLY,
                        region_bedfile=regions,
                        max_amplicons_region_group=limit,
                    )
                    n = len(amplicons(out))
                    self.assertGreater(n, 0)
                    self.assertLessEqual(n, expected_max)
                    if limit is None:
                        unlimited = n
            self.assertGreater(unlimited, 2)


class TestReplaceModes(unittest.TestCase):
    def config(self) -> Config:
        config = Config(**json.loads((SCHEME_DIR / "config.json").read_text()))
        config.amplicon_size_max = 900
        return config

    def run_replace(self, out: pathlib.Path, primername: str, **kwargs):
        replace(
            config=self.config(),
            primerbed=SCHEME_DIR / "primer.bed",
            primername=primername,
            msapath=SCHEME_DIR / "work" / "reference.fasta",
            pm=None,
            output=out,
            force=True,
            **kwargs,
        )

    def test_listall(self):
        with tempfile.TemporaryDirectory(dir="tests/integration") as d:
            out = pathlib.Path(d)
            self.run_replace(out, "1540b43d_11_RIGHT_1", mode=ReplaceRunModes.ListAll)
            options = [
                line
                for line in (out / "primer.bed").read_text().splitlines()
                if line.startswith("# PrimerPair Option:")
            ]
            self.assertGreater(len(options), 1)

    def test_no_mask_old_sites(self):
        with tempfile.TemporaryDirectory(dir="tests/integration") as d:
            out = pathlib.Path(d)
            self.run_replace(
                out,
                "1540b43d_11_RIGHT_1",
                mode=ReplaceRunModes.AddBest,
                mask_old_sites=False,
            )
            self.assertTrue((out / "primer.bed").stat().st_size > 0)

    def test_chrom_mismatch_is_a_usage_error(self):
        with tempfile.TemporaryDirectory(dir="tests/integration") as d:
            d = pathlib.Path(d)
            bed = d / "primer.bed"
            bed.write_text(
                (SCHEME_DIR / "primer.bed")
                .read_text()
                .replace("MN908947.3", "other_chrom")
            )
            with self.assertRaises(UsageError):
                replace(
                    config=self.config(),
                    primerbed=bed,
                    primername="1540b43d_11_RIGHT_1",
                    msapath=SCHEME_DIR / "work" / "reference.fasta",
                    pm=None,
                    output=d / "out",
                    force=True,
                    mode=ReplaceRunModes.AddBest,
                )


class TestVisualisePrimerDct(unittest.TestCase):
    def test_heatmap(self):
        array, seqdict = parse_msa(SCHEME_DIR / "reference.fasta")
        for include_seqs in (True, False):
            with self.subTest(include_seqs=include_seqs):
                html = primer_dct_heatmap(
                    array=array,
                    seqdict=seqdict,
                    bedfile=SCHEME_DIR / "primer.bed",
                    config=Config(),
                    include_seqs=include_seqs,
                    offline_plots=False,
                    mapping=MappingType.FIRST,
                )
                self.assertIn("plotly", html.lower())


if __name__ == "__main__":
    unittest.main()
