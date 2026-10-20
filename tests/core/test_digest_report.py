import pathlib
import tempfile
import unittest

from rich.console import Console

from primalscheme3.core import downsample as ds
from primalscheme3.core.config import Config
from primalscheme3.core.digest_report import failure_table, site_failures
from primalscheme3.core.logger import setup_rich_logger
from primalscheme3.core.msa import MSA, digest_kwargs
from primalscheme3.core.progress_tracker import ProgressManager
from tests.core.test_downsample import (
    BASE,
    FakeReturnAllKmer,
    FakeScorerAlwaysFail,
    FakeScorerAlwaysPass,
    variant,
)

MSA_PATH = pathlib.Path("tests/core/test_mismatch.fasta")


class ReturnAllKmer(FakeReturnAllKmer):
    """As FakeReturnAllKmer, with passed() as the Rust kmer computes it."""

    def passed(self):
        kept = [t for t in self._statuses if "LowFreq" not in t]
        return (
            bool(kept)
            and all(not t for t in kept)
            and not self._dimers
            and not any(
                "LowFreq" not in t and "Ignored" not in t for t, _ in self._errors
            )
        )


def render(table) -> str:
    console = Console(width=100, record=True, file=open("/dev/null", "w"))
    console.print(table)
    # Without the cells' trailing padding
    return "\n".join(line.rstrip() for line in console.export_text().splitlines())


def checks(*names):
    return frozenset(names)


def kmer_key(kmers, pos):
    return {
        (pos(k), tuple(sorted(zip(k.seqs(), k.counts(), strict=True)))) for k in kmers
    }


class TestSiteFailures(unittest.TestCase):
    def test_passing_site(self):
        self.assertEqual(site_failures(ReturnAllKmer([BASE], [1.0], [[]])), checks())

    def test_every_check_failed(self):
        kmer = ReturnAllKmer(
            [BASE, variant(3), variant(5)],
            [10.0, 3.0, 2.0],
            [[], ["Hairpin"], ["LowTm", "Hairpin"]],
        )
        self.assertEqual(site_failures(kmer), checks("Hairpin", "LowTm"))

    def test_below_min_freq_is_not_counted(self):
        kmer = ReturnAllKmer(
            [BASE, variant(3)],
            [10.0, 5.0],
            [["LowGC"], ["Hairpin", "LowFreq"]],
            errors=[
                (["ContainsN", "LowFreq"], 8.0),
                (["GapOnSetBase", "Ignored"], 9.0),
            ],
        )
        self.assertEqual(site_failures(kmer), checks("LowGC"))

    def test_errors_and_dimers(self):
        kmer = ReturnAllKmer(
            [BASE, variant(3)],
            [4.0, 1.0],
            [["HighGC"], []],
            errors=[(["ContainsN"], 3.0)],
            dimers=[(0, 1)],
        )
        self.assertEqual(
            site_failures(kmer), checks("HighGC", "ContainsN", "PrimerDimer")
        )

    def test_only_low_freq(self):
        kmer = ReturnAllKmer([BASE], [1.0], [["LowFreq"]])
        self.assertEqual(site_failures(kmer), checks("OnlyLowFreq"))


class TestFailureTable(unittest.TestCase):
    def test_fails_at_and_only_failure(self):
        rejected = (
            [checks("HighGC")] * 50
            + [checks("HighGC", "HighAnnealing")] * 10
            + [checks("PrimerDimer")] * 30
            + [checks("ContainsN", "HighGC")] * 9
            + [checks("ToLong")]
        )
        text = render(failure_table("chrom", rejected, 400, caption="a caption"))
        self.assertIn("chrom: 300 of 400 primer sites usable (75%)", text)
        self.assertIn("a caption", text)
        # Fails at, %, and Only failure, independent of check order
        self.assertRegex(text, r"HighGC\s+69\s+69%.*\s50\n")
        self.assertRegex(text, r"HighAnnealing\s+10\s+10%.*\s0\n")
        self.assertRegex(text, r"PrimerDimer\s+30\s+30%.*\s30\n")
        self.assertRegex(text, r"ContainsN\s+9\s+9%.*\s0\n")
        self.assertRegex(text, r"ToLong\s+1\s+1%.*\s1\n")
        self.assertRegex(text, r"Rejected sites\s+100\s+81\n")
        # Thermodynamic, interaction and alignment checks in that order
        self.assertLess(text.index("HighAnnealing"), text.index("PrimerDimer"))
        self.assertLess(text.index("PrimerDimer"), text.index("ContainsN"))

    def test_nothing_rejected(self):
        text = render(failure_table("chrom", [], 10))
        self.assertIn("10 of 10 primer sites usable (100%)", text)

    def test_logged_as_a_table(self):
        with tempfile.TemporaryDirectory() as d:
            logfile = pathlib.Path(d) / "file.log"
            logger = setup_rich_logger(str(logfile))
            logger.info(failure_table("chrom", [checks("Hairpin")] * 2, 5))
            for handler in logger.handlers:
                handler.flush()
            text = logfile.read_text()
        self.assertIn("3 of 5 primer sites usable", text)
        self.assertRegex(text, r"Hairpin\s+2\s+100%")
        self.assertNotIn("rich.table", text)


class TestDigestReport(unittest.TestCase):
    def digest(self, **cfg):
        config = Config(ncores=1, **cfg)
        msa = MSA("test", MSA_PATH, 0, ProgressManager(), config, logger=None)
        msa.digest_rs(config=config)
        return msa, config

    def test_kmers_match_a_normal_digest(self):
        # Digesting with return_all keeps exactly what a normal digest keeps
        for min_base_freq in (0.0, 0.1):
            with self.subTest(min_base_freq=min_base_freq):
                msa, config = self.digest(min_base_freq=min_base_freq)
                f, r, _ = msa._digester.digest(
                    findexes=None,
                    rindexes=None,
                    thermo_check=True,
                    return_all=False,
                    **digest_kwargs(config),
                )

                self.assertEqual(
                    kmer_key(msa.fkmers, lambda k: k.end), kmer_key(f, lambda k: k.end)
                )
                self.assertEqual(
                    kmer_key(msa.rkmers, lambda k: k.start),
                    kmer_key(r, lambda k: k.start),
                )
                for k in msa.fkmers + msa.rkmers:
                    self.assertTrue(k.passed())

    def test_report_accounts_for_every_site(self):
        import logging

        logger = logging.getLogger("test_digest_report")
        for downsample in (False, True):
            config = Config(ncores=1, downsample=downsample)
            msa = MSA("test", MSA_PATH, 0, ProgressManager(), config, logger=logger)
            with self.subTest(downsample=downsample), self.assertLogs(logger) as logs:
                msa.digest_rs(config=config)
            table = next(r.msg for r in logs.records if not isinstance(r.msg, str))
            text = render(table)
            usable, sites = (
                int(x.replace(",", ""))
                for x in text.split("primer sites usable")[0]
                .split(":")[-1]
                .split(" of ")
            )
            self.assertEqual(usable, len(msa.fkmers) + len(msa.rkmers))
            rejected = int(text.split("Rejected sites")[-1].split()[0].replace(",", ""))
            self.assertEqual(usable + rejected, sites)


class TestDownsampleFailure(unittest.TestCase):
    """A cloud downsampling can't use names the checks involved."""

    def setUp(self):
        self.config = Config(
            downsample=True, downsample_target=0.95, downsample_always_add_prop=0.5
        )

    def run_ds(self, seqs, scorer, fail_tags=None, dimers=None):
        return ds.delta_ct_downsample(
            seqs,
            self.config,
            scorer,
            failed_seqs=set(fail_tags or {}),
            dimers=dimers,
            fail_tags=fail_tags,
        )

    def test_cloud_failures(self):
        ok = FakeReturnAllKmer([BASE], [1.0], [[]])
        self.assertEqual(ds.cloud_failures(ok), checks())
        n = FakeReturnAllKmer(
            [BASE], [1.0], [[]], errors=[(["ContainsN"], 1.0), (["ToLong"], 1.0)]
        )
        self.assertEqual(ds.cloud_failures(n), checks("ContainsN", "ToLong"))
        low = FakeReturnAllKmer([BASE], [1.0], [["LowFreq"]])
        self.assertEqual(ds.cloud_failures(low), checks("OnlyLowFreq"))

    def test_kept(self):
        result = self.run_ds({BASE: 10.0}, FakeScorerAlwaysPass())
        self.assertEqual((result.kept, result.failures), ([BASE], checks()))

    def test_no_passing_primer(self):
        tags = {BASE: ["LowTm"], variant(3): ["Hairpin"], variant(5): ["Hairpin"]}
        result = self.run_ds(
            {BASE: 10.0, variant(3): 2.0, variant(5): 2.0}, FakeScorerAlwaysPass(), tags
        )
        self.assertEqual(
            (result.kept, result.failures), (None, checks("LowTm", "Hairpin"))
        )

    def test_below_target_names_the_uncovered(self):
        # BASE passes and is kept; only the uncovered variant's checks count
        result = self.run_ds(
            {BASE: 10.0, variant(3): 2.0},
            FakeScorerAlwaysFail(),
            {variant(3): ["Hairpin", "LowGC"]},
        )
        self.assertEqual(
            (result.kept, result.failures), (None, checks("Hairpin", "LowGC"))
        )

    def test_below_target_from_dimer_avoidance(self):
        # variant(3) passes but dimerises with the always-kept BASE
        result = self.run_ds(
            {BASE: 10.0, variant(3): 2.0},
            FakeScorerAlwaysFail(),
            dimers=[(BASE, variant(3))],
        )
        self.assertEqual(result.failures, checks("PrimerDimer"))

    def test_kept_primers_dimerise(self):
        self.config.downsample_always_add_prop = 0.1
        result = self.run_ds(
            {BASE: 10.0, variant(3): 5.0},
            FakeScorerAlwaysPass(),
            dimers=[(BASE, variant(3))],
        )
        self.assertEqual(result.failures, checks("PrimerDimer"))

    def test_self_dimer_is_a_primer_dimer(self):
        kmer = FakeReturnAllKmer(
            [BASE, variant(3)], [5.0, 1.0], [[], []], dimers=[(1, 1)]
        )
        cloud = ds.downsample_cloud(kmer)
        self.assertEqual(cloud.fail_tags, {variant(3): ["PrimerDimer"]})

    def test_python_thermo_check_names_the_check(self):
        # Without the Rust statuses (downsample-existing-scheme)
        result = ds.delta_ct_downsample(
            {"A" * 22: 1.0}, self.config, FakeScorerAlwaysPass()
        )
        self.assertEqual(result.failures, checks("LOW_GC"))


if __name__ == "__main__":
    unittest.main()
