import json
import logging
import pathlib
import shutil
import tempfile
import unittest

from primalscheme3.core.config import Config
from primalscheme3.repair.repair import (
    primers_equal,
    repair,
    report_check,
    required_n_pools,
)

TEST_SCHEME = pathlib.Path("tests/test_data/test_scheme").absolute()
SCHEME_CONFIG = TEST_SCHEME / "config.json"
SCHEME_BED = TEST_SCHEME / "primer.bed"
REFERENCE = TEST_SCHEME / "work" / "reference.fasta"
CHROM = "MN908947.3"

# 1540b43d_1_LEFT_1 covers [25, 51); a SNP 3 bases from its 3' end
LEFT_1 = "AACAAACCAACCAACTTTCGATCTCT"
SNP_POS = 48

# Primers that pass thermo_check under the default and test-scheme configs
PASSING = "AGCGTGGTTATTGGATGGGTTTG"
# Known dimer pair (SARS-CoV-2 v3 18_LEFT / 76_RIGHT), both thermo-passing
DIMER_A = "TGGAAATACCCACAAGTTAATGGTTTAAC"
DIMER_B = "ACACCTGTGCCTGTTAAACCAT"


class _PP:
    """Minimal stand-in for a PrimerPair: only the 0-indexed ``pool`` matters."""

    def __init__(self, pool: int):
        self.pool = pool


class Test_RequiredNPools(unittest.TestCase):
    def test_contiguous_pools(self):
        # pools {0, 1} -> need 2 buckets
        pps = [_PP(0), _PP(1), _PP(0)]
        self.assertEqual(required_n_pools(pps, min_pools=2), 2)

    def test_sparse_pools_do_not_undersize(self):
        # Regression: pools {0, 2} (pool 1 empty) must yield 3 buckets so that
        # indexing pools[pp.pool] for pp.pool == 2 does not raise IndexError.
        pps = [_PP(0), _PP(2)]
        self.assertEqual(required_n_pools(pps, min_pools=2), 3)
        # The highest index is addressable in a list of this size.
        buckets: list[list] = [[] for _ in range(required_n_pools(pps, 2))]
        for pp in pps:
            buckets[pp.pool].append(pp)  # must not raise

    def test_min_pools_floor(self):
        # More configured pools than are used -> keep the configured floor.
        pps = [_PP(0)]
        self.assertEqual(required_n_pools(pps, min_pools=4), 4)

    def test_empty_primerpairs(self):
        self.assertEqual(required_n_pools([], min_pools=2), 2)


class TestPrimersEqual(unittest.TestCase):
    def test_identical(self):
        self.assertTrue(primers_equal(LEFT_1, LEFT_1))

    def test_substitution(self):
        self.assertFalse(primers_equal(LEFT_1, LEFT_1[:-1] + "A"))

    @unittest.expectedFailure
    def test_prefix_is_not_equal(self):
        # Bug: zip(strict=False) stops at the shorter sequence, so a primer
        # that is a prefix of another counts as already in the scheme
        self.assertFalse(primers_equal(LEFT_1[:-3], LEFT_1))


class TestReportCheck(unittest.TestCase):
    def setUp(self):
        self.config = Config()
        self.logger = logging.getLogger("test_repair.report_check")

    def check(self, seq, count=1.0, current=(), pool=(), min_base_freq=0.0):
        self.config.min_base_freq = min_base_freq
        with self.assertLogs(self.logger, level="INFO") as logs:
            ok = report_check(
                seq=seq,
                count=count,
                current_primer_seqs=set(current),
                seqs_bytes_in_pools=[[s.encode() for s in pool]],
                pool=0,
                dimerscore=self.config.dimer_score,
                logger=self.logger,
                config=self.config,
            )
        return ok, "\n".join(logs.output)

    def test_valid(self):
        ok, log = self.check(PASSING)
        self.assertTrue(ok)
        self.assertIn("Can be added", log)

    def test_thermo_failure_reports_reason(self):
        ok, log = self.check("A" * 22)
        self.assertFalse(ok)
        self.assertIn("LOW_GC", log)

    def test_already_present(self):
        ok, log = self.check(PASSING, current=[PASSING])
        self.assertFalse(ok)
        self.assertIn("In scheme", log)

    def test_minor_allele(self):
        ok, log = self.check(PASSING, count=0.05, min_base_freq=0.1)
        self.assertFalse(ok)
        self.assertIn("Minor allele", log)

    def test_interaction_with_pool(self):
        ok, log = self.check(DIMER_A, pool=[DIMER_B])
        self.assertFalse(ok)
        self.assertIn("Interaction with pool", log)

    def test_no_interaction_with_unrelated_pool(self):
        ok, _ = self.check(DIMER_A, pool=[PASSING])
        self.assertTrue(ok)

    def test_checks_run_in_order(self):
        # A thermo failure is reported even if the seq is also present
        ok, log = self.check("A" * 22, current=["A" * 22])
        self.assertFalse(ok)
        self.assertIn("LOW_GC", log)
        self.assertNotIn("In scheme", log)


def read_reference() -> str:
    return "".join(
        line.strip() for line in REFERENCE.read_text().splitlines() if line[0] != ">"
    )


def bed_seqs(bed: str) -> set[str]:
    return {line.split("\t")[6] for line in bed.splitlines() if line and line[0] != "#"}


def snp(seq: str, pos: int) -> str:
    alt = {"A": "G", "G": "A", "C": "T", "T": "C"}[seq[pos]]
    return seq[:pos] + alt + seq[pos + 1 :]


class _Capture(logging.Handler):
    def __init__(self):
        super().__init__(logging.DEBUG)
        self.lines: list[str] = []

    def emit(self, record):
        self.lines.append(record.getMessage())


def run_repair(
    msa_seqs: dict[str, str], outdir: pathlib.Path
) -> tuple[pathlib.Path, str]:
    """Run repair on an MSA of ``msa_seqs``; returns the output dir and its log.

    The log is captured from the root logger, as work/file.log is only
    written by the first repair in a process (see test_log_file_per_run).
    """
    msa = outdir / "msa.fasta"
    msa.write_text("".join(f">{k}\n{v}\n" for k, v in msa_seqs.items()))
    out = outdir / "out"
    capture = _Capture()
    root = logging.getLogger()
    level = root.level
    root.addHandler(capture)
    root.setLevel(logging.DEBUG)
    try:
        repair(SCHEME_CONFIG, msa, SCHEME_BED, out, force=True, pm=None)
    finally:
        root.removeHandler(capture)
        root.setLevel(level)
    return out, "\n".join(capture.lines)


class TestRepairIntegration(unittest.TestCase):
    """Repair the test scheme against an MSA with a variant in 1_LEFT."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = pathlib.Path(tempfile.mkdtemp(dir="tests/repair"))
        ref = read_reference()
        cls.variant = snp(ref, SNP_POS)
        # The variant primer the digest walks back from the 1_LEFT end
        cls.variant_primer_3p = cls.variant[SNP_POS - 10 : 51]
        cls.out, cls.log = run_repair(
            {CHROM: ref, "var1": cls.variant, "var2": cls.variant}, cls.tmp
        )
        cls.bed = (cls.out / "primer.bed").read_text()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_outputs_written(self):
        for name in ("primer.bed", "config.json", "work/file.log", "work/msa.fasta"):
            self.assertTrue((self.out / name).is_file(), name)

    def test_every_primer_checked(self):
        n_lines = sum(
            1
            for line in SCHEME_BED.read_text().splitlines()
            if line[:1] not in ("#", "")
        )
        self.assertEqual(self.log.count("Checking "), n_lines)

    def test_existing_primers_reported_present(self):
        self.assertIn(LEFT_1, self.log)
        self.assertIn("In scheme", self.log)

    def test_variant_reported_valid(self):
        self.assertIn(self.variant_primer_3p, self.log, "variant primer not reported")
        self.assertIn(
            "Can be added", self.log.split(self.variant_primer_3p, 1)[1][:300]
        )

    def test_original_primers_kept_with_new_prefix(self):
        self.assertTrue(bed_seqs(SCHEME_BED.read_text()) <= bed_seqs(self.bed))
        self.assertNotIn("1540b43d_1_LEFT", self.bed)

    @unittest.expectedFailure
    def test_valid_variant_added_to_bed(self):
        # Bug: valid seqs are only added to the pool list, never to the
        # primerpair, so primer.bed never gains them
        self.assertIn(self.variant_primer_3p, self.bed)

    @unittest.expectedFailure
    def test_config_keeps_msa_data(self):
        # Bug: config.json is written from config.to_dict(), dropping the
        # msa_data entry repair adds for the new MSA
        cfg = json.loads((self.out / "config.json").read_text())
        self.assertEqual(len(cfg.get("msa_data") or {}), 2)


class TestRepairEdgeCases(unittest.TestCase):
    def setUp(self):
        self.tmp = pathlib.Path(tempfile.mkdtemp(dir="tests/repair"))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_rare_variant_not_offered(self):
        # Below min_base_freq (repair uses 0.1) the digest drops the variant
        ref = read_reference()
        variant = snp(ref, SNP_POS)
        seqs = {CHROM: ref, **{f"ref{i}": ref for i in range(10)}, "var": variant}
        _, log = run_repair(seqs, self.tmp)
        self.assertNotIn(variant[SNP_POS - 10 : 51], log)
        self.assertNotIn("Can be added", log)

    def test_identical_genomes_change_nothing(self):
        ref = read_reference()
        _, log = run_repair({CHROM: ref, "copy": ref}, self.tmp)
        self.assertIn("In scheme", log)
        self.assertNotIn("Can be added", log)
        self.assertNotIn("failed", log)

    def test_missing_reference_chrom_exits(self):
        ref = read_reference()
        with self.assertRaises(SystemExit):
            run_repair({"not_the_ref": ref, "other": ref}, self.tmp)

    @unittest.expectedFailure
    def test_log_file_per_run(self):
        # Bug: setup_rich_logger uses logging.basicConfig, which only takes
        # effect once per process, so later runs leave work/file.log empty
        ref = read_reference()
        for name in ("first", "second"):
            outdir = self.tmp / name
            outdir.mkdir()
            out, _ = run_repair({CHROM: ref}, outdir)
            self.assertIn("Checking", (out / "work" / "file.log").read_text(), name)

    def test_existing_output_needs_force(self):
        from click import UsageError

        ref = read_reference()
        out, _ = run_repair({CHROM: ref}, self.tmp)
        with self.assertRaises(UsageError):
            repair(
                SCHEME_CONFIG,
                self.tmp / "msa.fasta",
                SCHEME_BED,
                out,
                force=False,
                pm=None,
            )


if __name__ == "__main__":
    unittest.main()
