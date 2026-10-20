"""
Bundled primer scorer models, which predict the ΔCt of a primer with one
mismatch to its template.

Each model is a pair of files: ``<id>.joblib`` (a scikit-learn
HistGradientBoostingRegressor) and ``<id>.json`` (its manifest). The manifest
holds the feature encoding, the conditions the thermodynamic features are
computed at, the model file's sha256, the training data's provenance and
fit, and golden predictions that the tests check. Models are made by
scripts/train_primer_scorer.py; see docs/primer-scorer-models.md.

Versioning: a released model file is never changed. Retraining gives a new
id (``q5_dct_v2``), and old ids stay bundled so old schemes reproduce.
"""

import hashlib
import json
import pathlib

MODELS_DIR = pathlib.Path(__file__).parent
DEFAULT_MODEL = "q5_dct_v1"


class ModelNotFound(ValueError):
    pass


def available_models() -> list[str]:
    """Ids of the bundled models."""
    return sorted(p.stem for p in MODELS_DIR.glob("*.json"))


def resolve_model(model: str | pathlib.Path) -> pathlib.Path:
    """
    The manifest path for a bundled model id, or for a path to a model's
    ``.json`` manifest or ``.joblib`` file (the other must sit next to it).
    """
    if str(model) in available_models():
        return MODELS_DIR / f"{model}.json"
    path = pathlib.Path(model)
    manifest = path.with_suffix(".json")
    if path.suffix in {".json", ".joblib"} and manifest.is_file():
        return manifest
    raise ModelNotFound(
        f"Unknown primer scorer model '{model}'. Bundled models: "
        f"{', '.join(available_models())}; or give the path to a model's "
        ".json manifest."
    )


def load_manifest(model: str | pathlib.Path) -> tuple[dict, pathlib.Path]:
    """The manifest of a model, and the path of its model file."""
    manifest_path = resolve_model(model)
    manifest = json.loads(manifest_path.read_text())
    return manifest, manifest_path.parent / manifest["model_file"]


def file_sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()
