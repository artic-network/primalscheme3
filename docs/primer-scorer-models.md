# Primer scorer models

Downsampling (`--downsample`) and `visualise-primer-mismatches` predict how much a single mismatch between a primer and its template delays amplification (ΔCt). The prediction comes from a primer scorer model bundled with PrimalScheme3, so no model file or download is needed.

| Model id | Polymerase | Notes |
|---|---|---|
| `q5_dct_v1` | Q5 | Default. An XGBoost surrogate of a GAM fitted to qPCR data. |

Choose a model with `--primer-scorer-model <id>`, or give the path to a model's `.json` manifest. The `PRIMALSCHEME3_PRIMER_SCORER_MODEL` environment variable overrides both. The model used and its sha256 are recorded in the scheme's `config.json` (`primer_scorer_model`, `primer_scorer_model_sha256`).

## Files

Each model is two files in `primalscheme3/models/`:

- `<id>.ubj`: the booster, in XGBoost's own format. It loads in later XGBoost versions and loading it runs no code.
- `<id>.json`: the manifest, containing:
  - the feature encoding;
  - the thermodynamic conditions the features are computed at;
  - the model file's sha256;
  - the training data's sha256 and the software versions;
  - golden predictions.

The features are always computed at the manifest's conditions (for `q5_dct_v1`: 100 mM monovalent, 2 mM Mg, 0.8 mM dNTP, 15 nM DNA, 65 °C). These are the conditions the model was trained at, whatever the run's settings.

## Keeping results reproducible

- A released model is never changed. A retrained model gets a new id (`q5_dct_v2`), and earlier ids stay bundled so older schemes reproduce.
- `tests/core/test_primer_scorer_model.py` checks every bundled model against its manifest's golden predictions. An XGBoost upgrade, a change to feature extraction, or a change to the thermodynamics library that alters a prediction fails CI.

## Training a new model

1. Generate the training CSV from the GAM in qpcr-single-mm. Archive it alongside the release.
2. `uv run --group train python scripts/train_primer_scorer.py data.csv q5_dct_v2.joblib`
3. `uv run --group train python scripts/export_primer_scorer.py --joblib q5_dct_v2.joblib --training-csv data.csv --id q5_dct_v2`

   The export refuses to overwrite an existing id. It fails unless the bundled model reproduces the trained pipeline on the training data.
4. To make it the default, update `DEFAULT_MODEL` in `primalscheme3/models/__init__.py`.
