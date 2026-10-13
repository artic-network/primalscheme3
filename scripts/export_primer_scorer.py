"""
Export a trained primer scorer pipeline (joblib, from train_primer_scorer.py)
to a bundled model: ``primalscheme3/models/<id>.ubj`` and its ``<id>.json``
manifest.

The manifest records the feature encoding, the conditions the thermodynamic
features were computed at for training, provenance, and golden predictions
made by the original pipeline. The export fails unless the bundled model
reproduces the pipeline on the training data.

    uv run --group train python scripts/export_primer_scorer.py \
        --joblib q5_gam_surrogate_training_data.joblib \
        --training-csv q5_gam_surrogate_training_data.csv \
        --id q5_dct_v1

Never re-export over a released id: give a retrained model a new id.
"""

import argparse
import hashlib
import json
import pathlib
import random

import joblib
import numpy as np
import pandas as pd
import sklearn
import xgboost as xgb

from primalscheme3.core.primerscorer import PrimerScorer, extract_features
from primalscheme3.models import MODELS_DIR

# The conditions the training features (dg, prdg, p3_dg) were computed at:
# primalscheme3's Config() concentrations at 65 °C, as in qpcr-single-mm's
# create_feature_csv.py.
THERMO = {
    "mv_conc": 100.0,
    "dv_conc": 2.0,
    "dntp_conc": 0.8,
    "dna_conc": 15.0,
    "temp_c": 65.0,
}


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def golden_pairs(n: int, seed: int = 0) -> list[tuple[str, str]]:
    """(template, primer) pairs differing by one substitution, spread along
    the primer, with the 3' end weighted as in real mismatches."""
    rng = random.Random(seed)
    pairs = []
    while len(pairs) < n:
        query = "".join(rng.choice("ACGT") for _ in range(rng.randint(18, 32)))
        idx = min(int(rng.expovariate(0.25)), len(query) - 1)
        pos = len(query) - 1 - idx
        sub = rng.choice([b for b in "ACGT" if b != query[pos]])
        ref = query[:pos] + sub + query[pos + 1 :]
        if extract_features(ref, query, THERMO) is not None:
            pairs.append((ref, query))
    return pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--joblib", type=pathlib.Path, required=True)
    parser.add_argument("--training-csv", type=pathlib.Path, required=True)
    parser.add_argument("--id", required=True)
    parser.add_argument("--polymerase", default="Q5")
    args = parser.parse_args()

    model_file = MODELS_DIR / f"{args.id}.ubj"
    manifest_file = MODELS_DIR / f"{args.id}.json"
    if model_file.exists() or manifest_file.exists():
        raise SystemExit(f"{args.id} already exists: released models are immutable")

    pipeline = joblib.load(args.joblib)
    pre = pipeline.named_steps["preprocessor"]
    regressor = pipeline.named_steps["regressor"]
    (_, _, numeric), (_, encoder, categorical) = pre.transformers_
    categories = {
        name: [c.item() if isinstance(c, np.generic) else c for c in cats]
        for name, cats in zip(categorical, encoder.categories_, strict=True)
    }

    regressor.get_booster().save_model(model_file)

    training = pd.read_csv(args.training_csv)
    pairs = golden_pairs(60)
    features = [extract_features(r, q, THERMO) for r, q in pairs]
    expected = pipeline.predict(pd.DataFrame(features))

    manifest = {
        "id": args.id,
        "format_version": 1,
        "description": "ΔCt of a primer with one mismatch to its template: "
        "an XGBoost surrogate of a GAM fitted to qPCR data",
        "polymerase": args.polymerase,
        "target": "delta_ct",
        "model_file": model_file.name,
        "model_sha256": sha256(model_file),
        "features": {"numeric": list(numeric), "categorical": categories},
        "thermo": THERMO,
        "training": {
            "data_file": args.training_csv.name,
            "data_sha256": sha256(args.training_csv),
            "rows": len(training),
            "trainer": "scripts/train_primer_scorer.py",
            "hyperparameters": {
                k: v
                for k, v in regressor.get_params().items()
                if k in {"n_estimators", "max_depth", "learning_rate", "random_state"}
            },
            "xgboost_version": xgb.__version__,
            "scikit_learn_version": sklearn.__version__,
            "source_joblib_sha256": sha256(args.joblib),
        },
        "golden": [
            {"ref": r, "query": q, "features": f, "delta_ct": float(d)}
            for (r, q), f, d in zip(pairs, features, expected, strict=True)
        ],
    }
    manifest_file.write_text(json.dumps(manifest, indent=1, ensure_ascii=False) + "\n")

    # The bundled model must reproduce the pipeline exactly
    scorer = PrimerScorer(manifest_file)
    sample = training.drop(columns=["predicted_delta_ct"]).sample(
        n=min(10_000, len(training)), random_state=0
    )
    ours = np.array(scorer.predict(sample.to_dict("records")))
    theirs = pipeline.predict(sample)
    worst = float(np.max(np.abs(ours - theirs)))
    golden = np.array(scorer.predict(features))
    worst_golden = float(np.max(np.abs(golden - expected)))
    if worst > 1e-5 or worst_golden > 1e-5:
        model_file.unlink()
        manifest_file.unlink()
        raise SystemExit(f"Bundled model differs from the pipeline by {worst}")
    print(f"Wrote {model_file} and {manifest_file}; max |Δ| vs pipeline {worst:.2e}")


if __name__ == "__main__":
    main()
