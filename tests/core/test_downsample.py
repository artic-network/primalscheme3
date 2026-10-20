import csv
import pathlib
import tempfile
import unittest

import networkx as nx
from primalschemers import FKmer, RKmer

from primalscheme3.core import downsample as ds
from primalscheme3.core.config import Config
from primalscheme3.core.msa import MSA
from primalscheme3.core.primerscorer import PrimerScorer
from primalscheme3.core.progress_tracker import ProgressManager
from primalscheme3.core.seq_functions import primer_edit_dist
from primalscheme3.models import DEFAULT_MODEL

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


class RecordingScorer:
    # Records the (ref, query) pairs it is asked to score; mismatches harmless
    def __init__(self):
        self.pairs = []

    def score_seq_batch(self, pairs, config):
        self.pairs.extend(pairs)
        return [0.0] * len(pairs)


# A thermo-passing primer and single-substitution variants of it
BASE = "GGTATAGTCAATTGCCTAGTCTTCTGTTC"


def variant(pos_from_3p: int, seq: str = BASE) -> str:
    i = len(seq) - 1 - pos_from_3p
    alt = {"A": "G", "G": "A", "C": "T", "T": "C"}[seq[i]]
    return seq[:i] + alt + seq[i + 1 :]


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


class TestUpdateDeltaCt(unittest.TestCase):
    def test_none_is_ignored(self):
        self.assertEqual(ds.update_delta_ct({"A": -1.0}, {"A": None}), {"A": -1.0})

    def test_new_and_unknown_keys_take_value(self):
        old = {"A": None}
        self.assertEqual(
            ds.update_delta_ct(old, {"A": -2.0, "B": -5.0}), {"A": -2.0, "B": -5.0}
        )

    def test_keeps_the_best_score(self):
        # The least-penalised (highest) dCt wins
        self.assertEqual(ds.update_delta_ct({"A": -4.0}, {"A": -1.0}), {"A": -1.0})
        self.assertEqual(ds.update_delta_ct({"A": -1.0}, {"A": -4.0}), {"A": -1.0})


class TestAnnealingHelpers(unittest.TestCase):
    def setUp(self):
        self.scores = {"A": {"B": 5.0}, "B": {"A": 5.0}}
        self.counts = {"A": 10.0, "B": 2.0, "C": 1.0}

    def test_update_annealing_adds_self_and_neighbours(self):
        ann = ds.update_annealing("A", {"A": 0.0, "B": 0.0, "C": 0.0}, self.scores)
        self.assertEqual(ann, {"A": ds.NAIVE_ANNEALING, "B": 5.0, "C": 0.0})

    def test_add_annealing_uses_graph_neighbours(self):
        g = nx.Graph()
        g.add_edge("A", "B")
        ann = ds.add_annealing(g, "A", self.scores, {"A": 0.0, "B": 0.0})
        self.assertEqual(ann, {"A": ds.NAIVE_ANNEALING, "B": 5.0})

    def test_calc_scores_excludes_failed_and_sorts(self):
        scores = ds.calc_scores(self.scores, self.counts, {}, failed_seqs={"C"})
        self.assertNotIn("C", scores)
        self.assertEqual(list(scores), ["A", "B"])
        # A: naive self annealing x10 + B's 5 x 2
        self.assertEqual(scores["A"], ds.NAIVE_ANNEALING * 10 + 5.0 * 2)

    def test_calc_scores_saturates(self):
        # Once fully annealed, a seq gains nothing from itself
        full = {"A": ds.NAIVE_ANNEALING, "B": ds.NAIVE_ANNEALING}
        scores = ds.calc_scores(self.scores, self.counts, full, failed_seqs=set())
        self.assertEqual(scores["A"], 0)


class TestDeltaCtDownsampleBehaviour(unittest.TestCase):
    def setUp(self):
        self.config = Config(
            downsample=True,
            downsample_always_add_prop=0.5,
            downsample_target=0.9,
            downsample_max_dct=-4.0,
        )

    def run_ds(self, seq_counts, ps, **kw):
        return ds.delta_ct_downsample_seqs(
            seq_counts=seq_counts, config=self.config, ps=ps, visualise=False, **kw
        )

    def test_only_edit_distance_one_pairs_are_scored(self):
        far = variant(2, variant(8))  # two substitutions from BASE
        seqs = {BASE: 10.0, variant(3): 1.0, far: 1.0}
        ps = RecordingScorer()
        self.run_ds(seqs, ps)
        self.assertTrue(ps.pairs)
        for ref, query in ps.pairs:
            self.assertEqual(primer_edit_dist(ref, query), 1)
        # Both directions of each edge
        self.assertIn((BASE, variant(3)), ps.pairs)
        self.assertIn((variant(3), BASE), ps.pairs)
        self.assertNotIn((BASE, far), ps.pairs)

    def test_unconnected_seqs_each_need_adding(self):
        # No edges: every sequence must be added to cover it
        a, b = BASE, variant(2, variant(8))
        result = self.run_ds({a: 5.0, b: 5.0}, RecordingScorer())
        self.assertEqual(result, sorted([a, b]))

    def test_unreachable_target_returns_none(self):
        # Every passing seq covers itself, so the target is only out of reach
        # when uncoverable (thermo-failing) seqs hold too many counts
        seqs = {"A" * 29: 3.0, BASE: 1.0, variant(3): 1.0, variant(5): 1.0}
        self.config.downsample_always_add_prop = 0.9
        self.assertIsNone(self.run_ds(seqs, FakeScorerAlwaysFail()))
        self.config.downsample_target = 0.5
        self.assertEqual(
            self.run_ds(seqs, FakeScorerAlwaysFail()),
            sorted([BASE, variant(3), variant(5)]),
        )

    def test_thermo_failures_are_never_returned(self):
        bad = "A" * 29
        # bad holds most of the counts and can't be added, so no cover
        self.assertIsNone(self.run_ds({bad: 100.0, BASE: 1.0}, FakeScorerAlwaysPass()))
        self.config.downsample_target = 0.0
        result = self.run_ds({bad: 100.0, BASE: 1.0}, FakeScorerAlwaysPass())
        self.assertNotIn(bad, result or [])

    def test_length_variants_are_connected(self):
        # A 5'-longer primer with one 3' substitution is still an edge
        longer = "AC" + variant(4)
        ps = RecordingScorer()
        self.run_ds({BASE: 10.0, longer: 1.0}, ps)
        self.assertIn((BASE, longer), ps.pairs)

    def test_tie_break_is_deterministic(self):
        seqs = {variant(3): 1.0, variant(5): 1.0, BASE: 1.0}
        self.config.downsample_always_add_prop = 1.0
        first = self.run_ds(dict(seqs), FakeScorerAlwaysPass())
        for order in (dict(reversed(list(seqs.items()))), dict(sorted(seqs.items()))):
            self.assertEqual(self.run_ds(order, FakeScorerAlwaysPass()), first)

    def test_summary_csv(self):
        with tempfile.TemporaryDirectory() as d:
            seqs = {BASE: 10.0, variant(3): 1.0}
            self.run_ds(seqs, FakeScorerAlwaysPass(), outdir=d, name="F_1")
            self.run_ds(seqs, FakeScorerAlwaysPass(), outdir=d, name="F_2")
            with open(pathlib.Path(d) / "downsample_summary.csv") as fh:
                rows = list(csv.DictReader(fh))
        self.assertEqual([r["name"] for r in rows], ["F_1", "F_2"])
        self.assertEqual(rows[0]["n_seqs_pre"], "2")
        self.assertEqual(rows[0]["n_seqs_post"], "1")
        self.assertEqual(rows[0]["passed"], "True")


class TestDownsampleWithBundledModel(unittest.TestCase):
    def setUp(self):
        self.config = Config(downsample=True, downsample_target=0.9)
        self.ps = PrimerScorer()

    def test_mild_3p_distal_mismatches_are_covered(self):
        # Mismatches far from the 3' end barely shift Ct, so the dominant
        # primer covers its variants
        seqs = {BASE: 20.0, variant(15): 1.0, variant(18): 1.0}
        result = ds.delta_ct_downsample_seqs(seqs, self.config, self.ps)
        self.assertEqual(result, [BASE])

    def test_kmer_wrapper_matches_seq_counts(self):
        seqs = [BASE, variant(1), variant(15)]
        counts = [20.0, 5.0, 1.0]
        fk = FKmer([s.encode() for s in seqs], 100, counts)
        rk = RKmer([s.encode() for s in seqs], 100, counts)
        direct = ds.delta_ct_downsample_seqs(
            dict(zip(seqs, counts, strict=True)), self.config, self.ps
        )
        self.assertEqual(ds.downsample_kmer(fk, self.config, self.ps), direct)
        self.assertEqual(ds.downsample_kmer(rk, self.config, self.ps), direct)

    def test_pool_worker_matches_direct_call(self):
        seqs = [BASE, variant(1), variant(15)]
        counts = [20.0, 5.0, 1.0]
        ds._ds_pool_init(DEFAULT_MODEL, self.config)
        idx, kept = ds._ds_pool_worker((7, seqs, counts))
        self.assertEqual(idx, 7)
        direct = ds.delta_ct_downsample_seqs(
            dict(zip(seqs, counts, strict=True)), self.config, self.ps
        )
        self.assertEqual(kept, direct)


class TestAnnealingDownsample(unittest.TestCase):
    """downsample_seqs: the annealing-based selection used by downsample_scheme."""

    def setUp(self):
        self.config = Config(downsample_always_add_prop=0.5, downsample_target=0.9)

    def test_single_seq(self):
        self.assertEqual(ds.downsample_seqs({BASE: 1.0}, self.config), [BASE])

    def test_all_thermo_failures(self):
        self.assertIsNone(ds.downsample_seqs({"A" * 29: 1.0}, self.config))

    def test_similar_variant_is_covered_by_dominant(self):
        result = ds.downsample_seqs({BASE: 20.0, variant(15): 1.0}, self.config)
        self.assertIn(BASE, result)


class TestMSADownsample(unittest.TestCase):
    """End to end through MSA.digest_rs, serially and with the worker pool."""

    msa_path = pathlib.Path("tests/core/test_mismatch.fasta")

    def digest(self, ncores):
        config = Config(downsample=True, ncores=ncores)
        msa = MSA("test", self.msa_path, 0, ProgressManager(), config)
        msa.digest_rs(config=config)
        f = sorted((k.end, tuple(k.seqs())) for k in msa.fkmers)
        r = sorted((k.start, tuple(k.seqs())) for k in msa.rkmers)
        return f, r

    def test_parallel_matches_serial(self):
        serial = self.digest(ncores=1)
        self.assertTrue(serial[0] and serial[1])
        self.assertEqual(self.digest(ncores=2), serial)

    def test_downsampled_kmers_are_subsets(self):
        config = Config(downsample=False)
        msa = MSA("test", self.msa_path, 0, ProgressManager(), config)
        # As the downsample path: no Rust thermo, no Python hairpin filter
        msa.digest_rs(config=config, rs_thermo=False, py_hairpin=False)
        full = {k.end: set(k.seqs()) for k in msa.fkmers}
        fks, _ = self.digest(ncores=1)
        for end, seqs in fks:
            self.assertTrue(set(seqs) <= full[end], end)


def write_bed(path: pathlib.Path, rows, headers=("# artic-bed-version v3.0",)):
    """
    Write a bed of ``rows``: (chrom, amplicon, class, index, pool, seq, count),
    with ``pc=count`` unless count is None. Coordinates follow the amplicon.
    """
    offsets = {"LEFT": 0, "PROBE": 200, "RIGHT": 400}
    lines = list(headers)
    for chrom, amp, cls, idx, pool, seq, count in rows:
        start = amp * 1000 + offsets[cls]
        strand = "-" if cls == "RIGHT" else "+"
        line = f"{chrom}\t{start}\t{start + len(seq)}\tx_{amp}_{cls}_{idx}\t{pool}\t{strand}\t{seq}"
        if count is not None:
            line += f"\tpc={count}"
        lines.append(line)
    path.write_text("\n".join(lines) + "\n")
    return path


def side_seqs(scheme, chrom="chr1", amp=1, cls="LEFT"):
    return sorted(
        bl.sequence
        for bl in scheme.bedlines
        if bl.chrom == chrom
        and bl.amplicon_number == amp
        and bl.primer_class_str == cls
    )


class TestDownsampleScheme(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = pathlib.Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def run_ds(self, rows, visualise=False, **cfg):
        bed = write_bed(self.tmp / "primer.bed", rows)
        config = Config(downsample=True, **cfg)
        config.output = self.tmp / "out"
        return ds.downsample_scheme(bed, config, visualise=visualise)

    def pair(self, left, right=None, chrom="chr1", amp=1):
        """Rows for one amplicon from {seq: count} per side."""
        right = left if right is None else right
        return [
            (chrom, amp, "LEFT", i, 1, s, c) for i, (s, c) in enumerate(left.items(), 1)
        ] + [
            (chrom, amp, "RIGHT", i, 1, s, c)
            for i, (s, c) in enumerate(right.items(), 1)
        ]

    def test_single_primer_per_side_unchanged(self):
        rows = self.pair({BASE: 3.0})
        scheme = self.run_ds(rows)
        original = ds.Scheme.from_file(str(self.tmp / "primer.bed"))
        self.assertEqual(scheme.to_str(), original.to_str())

    def test_dominant_covers_distal_variants(self):
        scheme = self.run_ds(
            self.pair({BASE: 20.0, variant(15): 1.0, variant(18): 1.0})
        )
        self.assertEqual(side_seqs(scheme, cls="LEFT"), [BASE])
        self.assertEqual(side_seqs(scheme, cls="RIGHT"), [BASE])

    def test_counts_drive_the_choice(self):
        scheme = self.run_ds(self.pair({BASE: 1.0, variant(15): 20.0}))
        self.assertEqual(side_seqs(scheme), [variant(15)])

    def test_missing_counts_default_to_one(self):
        # Equal counts of 1: both clear always_add, so both are kept
        scheme = self.run_ds(self.pair({BASE: None, variant(15): None}))
        self.assertEqual(side_seqs(scheme), sorted([BASE, variant(15)]))

    def test_3p_terminal_variant_kept(self):
        # The dominant primer can't cover a 3'-terminal mismatch (ΔCt < -2)
        scheme = self.run_ds(self.pair({BASE: 20.0, variant(0): 1.0}))
        self.assertEqual(side_seqs(scheme), sorted([BASE, variant(0)]))

    def test_amplicons_and_chroms_are_independent(self):
        clouds = {
            ("chr1", 1): {BASE: 20.0, variant(15): 1.0},
            ("chr1", 2): {BASE: 20.0, variant(0): 1.0},
            ("chr2", 1): {BASE: 1.0, variant(15): 20.0},
        }
        rows = [
            r
            for (chrom, amp), c in clouds.items()
            for r in self.pair(c, chrom=chrom, amp=amp)
        ]
        together = self.run_ds(rows)
        for (chrom, amp), cloud in clouds.items():
            alone = self.run_ds(self.pair(cloud, chrom=chrom, amp=amp))
            self.assertEqual(
                side_seqs(together, chrom, amp), side_seqs(alone, chrom, amp)
            )

    def test_thermo_failures_never_kept(self):
        # A failing primer can't be covered, so allow for its 1/21 of counts
        scheme = self.run_ds(
            self.pair({BASE: 20.0, "A" * 29: 1.0}), downsample_target=0.9
        )
        self.assertEqual(side_seqs(scheme), [BASE])

    def test_unreachable_side_is_dropped(self):
        # Current behaviour: a cloud that can't reach the target loses that
        # side entirely, while the other side is kept
        scheme = self.run_ds(self.pair({"A" * 29: 100.0, BASE: 1.0}, {BASE: 1.0}))
        self.assertEqual(side_seqs(scheme, cls="LEFT"), [])
        self.assertEqual(side_seqs(scheme, cls="RIGHT"), [BASE])

    def test_uses_the_dct_scorer(self):
        from unittest import mock

        rec = RecordingScorer()
        with mock.patch.object(ds, "PrimerScorer", lambda model: rec):
            self.run_ds(
                self.pair({BASE: 20.0, variant(15): 1.0, variant(2, variant(8)): 1.0})
            )
        self.assertTrue(rec.pairs)
        for ref, query in rec.pairs:
            self.assertEqual(primer_edit_dist(ref, query), 1)
        # Unknown scores: each sequence has to cover itself
        with mock.patch.object(
            ds, "PrimerScorer", lambda model: FakeScorerAlwaysFail()
        ):
            scheme = self.run_ds(self.pair({BASE: 20.0, variant(15): 1.0}))
        self.assertEqual(side_seqs(scheme), sorted([BASE, variant(15)]))

    def test_model_choice(self):
        import os
        from unittest import mock

        from primalscheme3.models import ModelNotFound

        with self.assertRaises(ModelNotFound):
            self.run_ds(self.pair({BASE: 1.0}), primer_scorer_model="nope")
        with mock.patch.dict(os.environ, {"PRIMALSCHEME3_PRIMER_SCORER_MODEL": "nope"}):
            with self.assertRaises(ModelNotFound):
                self.run_ds(self.pair({BASE: 1.0}))

    def test_headers_kept(self):
        scheme = self.run_ds(self.pair({BASE: 1.0}))
        self.assertEqual(scheme.headers, ["# artic-bed-version v3.0"])

    def test_probes_kept(self):
        probe = "ACGTGCATGCAGTCAGTCAGTA"
        rows = self.pair({BASE: 20.0, variant(15): 1.0})
        rows.append(("chr1", 1, "PROBE", 1, 1, probe, None))
        scheme = self.run_ds(rows)
        self.assertEqual(side_seqs(scheme, cls="PROBE"), [probe])

    def test_one_sided_amplicon(self):
        rows = [
            ("chr1", 1, "LEFT", 1, 1, BASE, 20.0),
            ("chr1", 1, "LEFT", 2, 1, variant(15), 1.0),
        ]
        scheme = self.run_ds(rows)
        self.assertEqual(side_seqs(scheme), [BASE])
        self.assertEqual(side_seqs(scheme, cls="RIGHT"), [])

    def test_output_order_and_round_trip(self):
        rows = self.pair({BASE: 20.0, variant(0): 1.0}, amp=1) + self.pair(
            {BASE: 5.0}, amp=2
        )
        scheme = self.run_ds(rows)
        names = [bl.primername for bl in scheme.bedlines]
        self.assertEqual(
            names,
            [
                "x_1_LEFT_1",
                "x_1_LEFT_2",
                "x_1_RIGHT_1",
                "x_1_RIGHT_2",
                "x_2_LEFT_1",
                "x_2_RIGHT_1",
            ],
        )
        reparsed = ds.Scheme.from_str(scheme.to_str())
        self.assertEqual(reparsed.to_str(), scheme.to_str())
        # Attributes such as pc survive
        self.assertEqual(scheme.bedlines[0].attributes["pc"], "20.0")

    def test_visualise_writes_plots_and_summary(self):
        import matplotlib

        matplotlib.use("Agg")
        self.run_ds(self.pair({BASE: 20.0, variant(15): 1.0}), visualise=True)
        out = self.tmp / "out"
        self.assertTrue((out / "chr1_amp000_LEFT.svg").is_file())
        self.assertTrue((out / "chr1_amp000_RIGHT.svg").is_file())
        with open(out / "downsample_summary.csv") as fh:
            names = [r["name"] for r in csv.DictReader(fh)]
        self.assertEqual(names, ["chr1_amp000_LEFT", "chr1_amp000_RIGHT"])

    def test_cli(self):
        from typer.testing import CliRunner

        from primalscheme3.cli import app

        bed = write_bed(
            self.tmp / "primer.bed", self.pair({BASE: 20.0, variant(15): 1.0})
        )
        result = CliRunner().invoke(app, ["downsample-existing-scheme", str(bed)])
        self.assertEqual(result.exit_code, 0, result.output)
        self.assertIn("x_1_LEFT_1", result.output)
        self.assertNotIn(variant(15), result.output)
        bad = CliRunner().invoke(
            app,
            ["downsample-existing-scheme", str(bed), "--primer-scorer-model", "nope"],
        )
        self.assertNotEqual(bad.exit_code, 0)
