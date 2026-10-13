"""
Train a primer scorer model and bundle it in primalscheme3/models.

The model is a scikit-learn HistGradientBoostingRegressor, a surrogate of the
GAM fitted to qPCR data in qpcr-single-mm. It is trained on that GAM's
predicted ΔCt over a grid of mismatch features
(q5_gam_surrogate_training_data.csv), so primalscheme3 can predict ΔCt
without R.

Writes ``<id>.joblib`` and its ``<id>.json`` manifest: the feature encoding,
the thermodynamic conditions the features are computed at, checksums, the
training provenance and fit, and golden predictions that the tests check.
The run fails unless the bundled model reproduces the trained one exactly.

    uv run --group train python scripts/train_primer_scorer.py \\
        --training-csv q5_gam_surrogate_training_data.csv --id q5_dct_v2

A released model is never retrained in place: give a new model a new id.
``--replace`` overwrites an id, and is only for ids not yet released.
"""

import argparse
import hashlib
import json
import platform
import random
import time

import joblib
import numpy as np
import pandas as pd
import sklearn
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import train_test_split

from primalscheme3.core.primerscorer import PrimerScorer, extract_features
from primalscheme3.models import MODELS_DIR

# The conditions the training features (dg, prdg, p3_dg) were computed at:
# primalscheme3's Config() concentrations at 65 °C, as in qpcr-single-mm's
# create_feature_csv.py. PrimerScorer computes features at these too.
THERMO = {
    "mv_conc": 100.0,
    "dv_conc": 2.0,
    "dntp_conc": 0.8,
    "dna_conc": 15.0,
    "temp_c": 65.0,
}
NUMERIC = ["p3_dg", "prdg", "idx_prop", "mm_dg2"]
CATEGORICAL = ["mut_base", "dinu2_gc_rich"]
TARGET = "predicted_delta_ct"
HYPERPARAMETERS = {
    "max_iter": 500,
    "learning_rate": 0.1,
    "max_leaf_nodes": 63,
    "early_stopping": False,
    "random_state": 42,
}
TEST_SIZE = 0.2
SPLIT_SEED = 42


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def encode(df: pd.DataFrame, categories: dict[str, list]) -> np.ndarray:
    """Numeric features, then each categorical feature as its category index,
    exactly as PrimerScorer._encode does at runtime."""
    columns = [df[name].to_numpy(np.float64) for name in NUMERIC]
    for name, cats in categories.items():
        index = {c: i for i, c in enumerate(cats)}
        columns.append(df[name].map(index).to_numpy(np.float64))
    return np.column_stack(columns)


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
    parser.add_argument("--training-csv", required=True)
    parser.add_argument("--id", required=True)
    parser.add_argument("--polymerase", default="Q5")
    parser.add_argument(
        "--replace", action="store_true", help="Overwrite an unreleased id"
    )
    args = parser.parse_args()

    model_file = MODELS_DIR / f"{args.id}.joblib"
    manifest_file = MODELS_DIR / f"{args.id}.json"
    if (model_file.exists() or manifest_file.exists()) and not args.replace:
        raise SystemExit(f"{args.id} exists: released models are immutable")

    df = pd.read_csv(args.training_csv)
    categories = {
        name: sorted(
            v.item() if isinstance(v, np.generic) else v for v in df[name].unique()
        )
        for name in CATEGORICAL
    }
    X = encode(df, categories)
    y = df[TARGET].to_numpy()
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=TEST_SIZE, random_state=SPLIT_SEED
    )

    n_numeric = len(NUMERIC)
    model = HistGradientBoostingRegressor(
        categorical_features=list(range(n_numeric, n_numeric + len(CATEGORICAL))),
        **HYPERPARAMETERS,
    )
    start = time.time()
    model.fit(X_train, y_train)
    print(f"Trained in {time.time() - start:.0f}s")

    pred = model.predict(X_test)
    metrics = {
        "r2": float(r2_score(y_test, pred)),
        "rmse": float(mean_squared_error(y_test, pred) ** 0.5),
        "max_abs_error": float(np.max(np.abs(pred - y_test))),
    }
    print(f"Held-out fit: R² {metrics['r2']:.5f}, RMSE {metrics['rmse']:.4f} ΔCt")

    joblib.dump(model, model_file, compress=3)

    pairs = golden_pairs(60)
    features = [extract_features(r, q, THERMO) for r, q in pairs]
    expected = model.predict(encode(pd.DataFrame(features), categories))

    manifest = {
        "id": args.id,
        "format_version": 2,
        "description": "ΔCt of a primer with one mismatch to its template: a "
        "scikit-learn HistGradientBoostingRegressor surrogate of a GAM fitted "
        "to qPCR data",
        "polymerase": args.polymerase,
        "target": "delta_ct",
        "model_file": model_file.name,
        "model_sha256": sha256(model_file),
        "features": {"numeric": NUMERIC, "categorical": categories},
        "thermo": THERMO,
        "training": {
            "data_file": pathlib_name(args.training_csv),
            "data_sha256": sha256(args.training_csv),
            "rows": len(df),
            "target": TARGET,
            "test_size": TEST_SIZE,
            "split_seed": SPLIT_SEED,
            "estimator": "sklearn.ensemble.HistGradientBoostingRegressor",
            "hyperparameters": HYPERPARAMETERS,
            "held_out_metrics": metrics,
            "trainer": "scripts/train_primer_scorer.py",
            "scikit_learn_version": sklearn.__version__,
            "numpy_version": np.__version__,
            "python_version": platform.python_version(),
        },
        "golden": [
            {"ref": r, "query": q, "features": f, "delta_ct": float(d)}
            for (r, q), f, d in zip(pairs, features, expected, strict=True)
        ],
    }
    manifest_file.write_text(json.dumps(manifest, indent=1, ensure_ascii=False) + "\n")

    # The bundled model, loaded and encoded as at runtime, must reproduce the
    # trained model exactly
    scorer = PrimerScorer(manifest_file)
    rows = df.sample(n=min(10_000, len(df)), random_state=0)
    ours = np.array(scorer.predict(rows.drop(columns=[TARGET]).to_dict("records")))
    theirs = model.predict(encode(rows, categories))
    golden = np.array(scorer.predict(features))
    worst = max(np.max(np.abs(ours - theirs)), np.max(np.abs(golden - expected)))
    if worst > 0:
        model_file.unlink()
        manifest_file.unlink()
        raise SystemExit(f"Bundled model differs from the trained model by {worst}")
    print(f"Wrote {model_file} and {manifest_file}")


def pathlib_name(path: str) -> str:
    return path.replace("\\", "/").rsplit("/", 1)[-1]


if __name__ == "__main__":
    main()
