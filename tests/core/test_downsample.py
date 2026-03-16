import unittest

from primalscheme3.core import downsample as ds
from primalscheme3.core.config import Config

# These tests purley focus in the greedy downsample selection logic.
# We assume:
#   - PrimerScore is providing correct values
#   - ThermoCheckers are accurate
#   - EditDist calculates are correct


# Fake scores for dependency injection
class FakeScorerAlwaysPass:
    # This will always return 0. i.e. mismatch has no effect.
    def score_seq_batch(self, pairs, config):
        return [0.0] * len(pairs)


class FakeScorerAlwaysFail:
    # This will always return None. i.e. mismatch has unpredictable / massive effect
    def score_seq_batch(self, pairs, config):
        return [None] * len(pairs)


class FakeScorerConsistent:
    # This will always return X
    def score_seq_batch(self, pairs, config):
        return [-3] * len(pairs)


class TestDownSample(unittest.TestCase):
    def setUp(self):
        self.config = Config(
            downsample=True,
            downsample_always_add_prop=0.2,
            downsample_target=0.9,
            downsample_max_dct=-4.0,
        )

    def test_delta_ct_downsample_single_primer_covers_all(self):
        # All the mismatches have 0 effect.
        # Should return the dominate seq
        seq_counts = {
            "GGTATAGTCAATTGCCTAGTCTTCTGTTC": 10.0,
            "GGTATAGTCAATTGCCTAGTCTTCTGTTG": 1.0,
            "GGTATAGTCAATTGCCTAGTCTTCTGTTT": 1.0,
        }
        result = ds.delta_ct_downsample_seqs(
            seq_counts=seq_counts,
            config=self.config,
            ps=FakeScorerAlwaysPass(),  # type: ignore
            visualise=False,
        )
        self.assertEqual(result, ["GGTATAGTCAATTGCCTAGTCTTCTGTTC"])

    def test_delta_ct_downsample_thermo_fail_recovery(self):
        # All the mismatches have 0 effect.
        # This should not add the dominate seq.
        # deterministic first seq should be added
        seq_counts = {
            "TGTTC": 10.0,  # thermo fail
            "GGTATAGTCAATTGCCTAGTCTTCTGTTG": 1.0,
            "GGTATAGTCAATTGCCTAGTCTTCTGTTT": 1.0,
        }
        result = ds.delta_ct_downsample_seqs(
            seq_counts=seq_counts,
            config=self.config,
            ps=FakeScorerAlwaysPass(),  # type: ignore
            visualise=False,
        )
        self.assertEqual(result, ["GGTATAGTCAATTGCCTAGTCTTCTGTTG"])

    def test_delta_ct_downsample_always_add_prop(self):
        # All the mismatches have 0 effect.
        # All seqs are above the ALWAYS_ADD_PROP threshold.
        # All should be added
        seq_counts = {
            "GGTATAGTCAATTGCCTAGTCTTCTGTTC": 10.0,
            "GGTATAGTCAATTGCCTAGTCTTCTGTTG": 10,
            "GGTATAGTCAATTGCCTAGTCTTCTGTTT": 10,
        }
        result = ds.delta_ct_downsample_seqs(
            seq_counts=seq_counts,
            config=self.config,
            ps=FakeScorerAlwaysPass(),  # type: ignore
            visualise=False,
        )
        self.assertEqual(
            result,
            [
                "GGTATAGTCAATTGCCTAGTCTTCTGTTC",
                "GGTATAGTCAATTGCCTAGTCTTCTGTTG",
                "GGTATAGTCAATTGCCTAGTCTTCTGTTT",
            ],
        )

    def test_delta_ct_downsample_all_primer_covers_all(self):
        # All the mismatches have massive dCts
        # Should return all seqs
        seq_counts = {
            "GGTATAGTCAATTGCCTAGTCTTCTGTTC": 1.0,
            "GGTATAGTCAATTGCCTAGTCTTCTGTTG": 1.0,
            "GGTATAGTCAATTGCCTAGTCTTCTGTTT": 1.0,
        }
        result = ds.delta_ct_downsample_seqs(
            seq_counts=seq_counts,
            config=self.config,
            ps=FakeScorerAlwaysFail(),  # type: ignore
            visualise=False,
        )
        self.assertEqual(
            result,
            [
                "GGTATAGTCAATTGCCTAGTCTTCTGTTC",
                "GGTATAGTCAATTGCCTAGTCTTCTGTTG",
                "GGTATAGTCAATTGCCTAGTCTTCTGTTT",
            ],
        )

    def test_delta_ct_downsample_downsample_fail(self):
        # All the mismatches have massive dCts
        # High count seq is failing and cannot be recovered.
        # Should return None
        seq_counts = {
            "TGTTC": 10.0,  # thermo fail + high count
            "GGTATAGTCAATTGCCTAGTCTTCTGTTG": 1.0,
            "GGTATAGTCAATTGCCTAGTCTTCTGTTT": 1.0,
        }
        result = ds.delta_ct_downsample_seqs(
            seq_counts=seq_counts,
            config=self.config,
            ps=FakeScorerAlwaysFail(),  # type: ignore
            visualise=False,
        )
        self.assertIsNone(
            result,
        )

    def test_delta_ct_downsample_dct_check(self):
        # All the mismatches have -3 dCts lt the dCt threshold
        # Should return All seqs
        self.config.downsample_max_dct = 0
        seq_counts = {
            "GGTATAGTCAATTGCCTAGTCTTCTGTTC": 10.0,
            "GGTATAGTCAATTGCCTAGTCTTCTGTTG": 10,
            "GGTATAGTCAATTGCCTAGTCTTCTGTTT": 10,
        }
        result = ds.delta_ct_downsample_seqs(
            seq_counts=seq_counts,
            config=self.config,
            ps=FakeScorerConsistent(),  # type: ignore | returns -3
            visualise=False,
        )
        self.assertEqual(
            result,
            [
                "GGTATAGTCAATTGCCTAGTCTTCTGTTC",
                "GGTATAGTCAATTGCCTAGTCTTCTGTTG",
                "GGTATAGTCAATTGCCTAGTCTTCTGTTT",
            ],
        )


class TestDownsampleCoverageHelpers(unittest.TestCase):
    def setUp(self):
        self.seq_counts = {
            "A": 10.0,
            "B": 5.0,
            "C": 2.0,
            "D": 1.0,
        }
        self.seq_delta_ct = {
            "A": 0.0,
            "B": -4.0,
            "C": -4.1,
            "D": None,
        }

    def test_accounted_count_respects_threshold_and_none(self):
        result = ds.accounted_count(
            self.seq_counts,
            self.seq_delta_ct,
            max_dct=-4.0,
        )
        # A (0.0) and B (-4.0) count, C below threshold and D None do not
        self.assertEqual(result, 15.0)

    def test_coverage_prop_uses_total_count(self):
        result = ds.coverage_prop(
            self.seq_counts,
            self.seq_delta_ct,
            max_dct=-4.0,
        )
        # (10 + 5) / (10 + 5 + 2 + 1) = 15 / 18
        self.assertAlmostEqual(result, 15.0 / 18.0)

    def test_coverage_prop_zero_total(self):
        result = ds.coverage_prop({}, {}, max_dct=-4.0)
        self.assertEqual(result, 0.0)
