"""
Parity tests for the bundled primer scorer models.

Each model's manifest carries golden predictions made when the model was
trained. These must keep matching through scikit-learn upgrades, changes to
feature extraction, and changes to the thermodynamics library.
"""

import json
import shutil
import tempfile
import unittest

from primalscheme3.core import downsample as ds
from primalscheme3.core.config import Config
from primalscheme3.core.primerscorer import PrimerScorer, extract_features
from primalscheme3.models import (
    DEFAULT_MODEL,
    MODELS_DIR,
    ModelNotFound,
    available_models,
    file_sha256,
    load_manifest,
)


class TestBundledModels(unittest.TestCase):
    def test_default_model_is_bundled(self):
        self.assertIn(DEFAULT_MODEL, available_models())

    def test_model_files_match_manifests(self):
        for model in available_models():
            with self.subTest(model=model):
                manifest, model_file = load_manifest(model)
                self.assertEqual(manifest["id"], model)
                self.assertEqual(file_sha256(model_file), manifest["model_sha256"])

    def test_golden_features(self):
        # Feature extraction (incl. ΔG from the thermo library) is unchanged
        for model in available_models():
            manifest, _ = load_manifest(model)
            for g in manifest["golden"]:
                with self.subTest(model=model, query=g["query"]):
                    f = extract_features(g["ref"], g["query"], manifest["thermo"])
                    self.assertIsNotNone(f)
                    for name, value in g["features"].items():
                        if isinstance(value, float):
                            self.assertAlmostEqual(f[name], value, places=6)
                        else:
                            self.assertEqual(f[name], value)

    def test_golden_predictions(self):
        # End to end: (template, primer) -> ΔCt, as the training pipeline gave
        for model in available_models():
            scorer = PrimerScorer(model)
            golden = scorer.manifest["golden"]
            preds = scorer.predict(
                [scorer._get_features(g["ref"], g["query"]) for g in golden]
            )
            for g, pred in zip(golden, preds, strict=True):
                with self.subTest(model=model, query=g["query"]):
                    self.assertAlmostEqual(pred, g["delta_ct"], places=6)

    def test_predictions_ignore_run_conditions(self):
        # Features use the model's training conditions, not the run's config
        g = PrimerScorer().manifest["golden"][5]
        pair = [(g["ref"], g["query"])]
        default = PrimerScorer().score_seq_batch(pair, Config())
        other = Config(mv_conc=20.0, primer_annealing_tempc=72)
        self.assertEqual(PrimerScorer().score_seq_batch(pair, other), default)

    def test_identical_and_unscorable_pairs(self):
        scorer = PrimerScorer()
        self.assertEqual(scorer.score_seq("ACGTACGT", "ACGTACGT"), 0.0)
        self.assertIsNone(scorer.score_seq("ACGTACGTAA", "TTGTACGTAC"))


class TestModelResolution(unittest.TestCase):
    def test_unknown_id(self):
        with self.assertRaises(ModelNotFound) as cm:
            PrimerScorer("q5_dct_v999")
        self.assertIn(DEFAULT_MODEL, str(cm.exception))

    def test_model_by_path(self):
        with tempfile.TemporaryDirectory() as d:
            for suffix in (".json", ".joblib"):
                shutil.copy(MODELS_DIR / f"{DEFAULT_MODEL}{suffix}", d)
            for name in (f"{DEFAULT_MODEL}.json", f"{DEFAULT_MODEL}.joblib"):
                scorer = PrimerScorer(f"{d}/{name}")
                self.assertEqual(scorer.id, DEFAULT_MODEL)

    def test_tampered_model_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            shutil.copy(MODELS_DIR / f"{DEFAULT_MODEL}.joblib", d)
            manifest = json.loads((MODELS_DIR / f"{DEFAULT_MODEL}.json").read_text())
            manifest["model_sha256"] = "0" * 64
            manifest_path = f"{d}/{DEFAULT_MODEL}.json"
            with open(manifest_path, "w") as f:
                json.dump(manifest, f)
            with self.assertRaises(ValueError):
                PrimerScorer(manifest_path)


class TestConfig(unittest.TestCase):
    def test_default_and_legacy_key(self):
        self.assertEqual(Config().primer_scorer_model, DEFAULT_MODEL)
        # config.json files written before models were bundled
        legacy = Config(primer_scorer_model_path="/models/q5.json")
        self.assertEqual(legacy.primer_scorer_model, "/models/q5.json")
        self.assertIn("primer_scorer_model", Config().to_dict())

    def test_downsample_runs_with_bundled_model(self):
        # No model path or environment variable needed
        seqs = {
            "GGTATAGTCAATTGCCTAGTCTTCTGTTC": 10.0,
            "GGTATAGTCAATTGCCTAGTCTTCTGTTG": 1.0,
            "GGTATAGTCAATTGCCTAGTCTTCTGTTT": 1.0,
        }
        config = Config(downsample=True, downsample_target=0.9)
        kept = ds.delta_ct_downsample_seqs(
            seq_counts=seqs, config=config, ps=PrimerScorer(), visualise=False
        )
        self.assertTrue(kept)
        self.assertIn("GGTATAGTCAATTGCCTAGTCTTCTGTTC", kept)


if __name__ == "__main__":
    unittest.main()
