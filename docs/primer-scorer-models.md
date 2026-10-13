# Primer scorer models

Downsampling (`--downsample`, experimental) and `visualise-primer-dct` predict how much a single mismatch between a primer and its template delays amplification (ΔCt). The prediction comes from a primer scorer model bundled with PrimalScheme3, so no model file or download is needed.

| Model id | Polymerase | Model | Held-out fit to the GAM |
|---|---|---|---|
| `q5_dct_v1` | Q5 | scikit-learn `HistGradientBoostingRegressor` | R² 0.99957, RMSE 0.032 ΔCt |

Choose a model with `--primer-scorer-model <id>`, or give the path to a model's `.json` manifest. The `PRIMALSCHEME3_PRIMER_SCORER_MODEL` environment variable overrides both. The model used and its sha256 are recorded in the scheme's `config.json` (`primer_scorer_model`, `primer_scorer_model_sha256`).

## How the model is made

1. **GAM, fitted in R** (the qpcr-single-mm analysis, `dinucleotide/`). A generalised additive model is fitted to qPCR measurements of Q5 amplification with single-mismatch primers. The measured ΔCt is modelled from the mismatch's features (below).
2. **Training grid.** The GAM predicts ΔCt over a dense grid of feature values. That gives `q5_gam_surrogate_training_data.csv`: 1,440,000 rows of the six features plus `predicted_delta_ct`.
3. **Surrogate, in Python** (`scripts/train_primer_scorer.py`). A gradient-boosted regressor learns the GAM's predictions, so PrimalScheme3 can score primers quickly and without R.
   - **Estimator:** scikit-learn `HistGradientBoostingRegressor`.
   - **Hyperparameters:** 500 iterations, learning rate 0.1, 63 leaves per tree, no early stopping, `random_state=42`.
   - **Data split:** trained on 80% of the grid; the other 20% (`random_state=42`) measures the fit recorded in the manifest.

**Features**, computed by `primalscheme3.core.primerscorer.extract_features` for a primer with exactly one mismatch to its template:

| Feature | Meaning |
|---|---|
| `p3_dg` | ΔG of the primer's 3′ segment beyond the mismatch |
| `prdg` | ΔG of the primer against its own sequence with the mismatch corrected |
| `idx_prop` | Position of the mismatch from the 3′ end, as a fraction of the primer length |
| `mm_dg2` | Nearest-neighbour ΔG of the mismatch's 3′ flanking context (0 for a 3′-terminal mismatch) |
| `mut_base` | The substitution, e.g. `A>C` (categorical) |
| `dinu2_gc_rich` | Whether the 3′ flanking dinucleotide is GC-rich (categorical) |

The ΔG features are computed with primer3 at 100 mM monovalent, 2 mM Mg²⁺, 0.8 mM dNTP, 15 nM DNA and 65 °C. These are the conditions the training data was generated at. `PrimerScorer` always uses the conditions in the model's manifest, whatever the run's settings, so its inputs match what the model learned.

## Files

Each model is two files in `primalscheme3/models/`:

- **`<id>.joblib`:** the fitted `HistGradientBoostingRegressor`. It's loaded only after its sha256 matches the manifest.
- **`<id>.json`:** the manifest, containing:
  - `features`: the feature order, and the category lists. Categorical features are passed as their index in these lists; an unseen category is passed as missing.
  - `thermo`: the conditions the ΔG features are computed at.
  - `model_sha256`: the checksum of the model file.
  - `training`: the training data's file name, sha256 and row count; the split; the estimator and hyperparameters; the held-out R², RMSE and maximum error; and the scikit-learn, numpy and Python versions.
  - `golden`: 60 primer–template pairs, with their features and the ΔCt the model gave when it was trained.

Runtime needs only scikit-learn, joblib and numpy. scikit-learn's wheels include their own OpenMP runtime, so nothing extra is needed on macOS. Prediction runs single-threaded, because scikit-learn's per-call thread start-up costs far more than scoring one primer cloud, and downsampling already parallelises across clouds.

## Keeping results reproducible

- **Never change a released model.** A retrained model gets a new id (`q5_dct_v2`), and earlier ids stay bundled so older schemes reproduce. The training script refuses to overwrite an existing id. `--replace` exists only for ids that haven't been released.
- **scikit-learn is pinned to the minor version the models were saved with** (`scikit-learn>=1.8,<1.9`). Model files are joblib (pickle), which scikit-learn only guarantees to load with the version that wrote them.
- **The golden tests** (`tests/core/test_primer_scorer_model.py`) check every bundled model's predictions and features against its manifest. A scikit-learn upgrade, a change to feature extraction, or a change to the thermodynamics library that alters a prediction fails CI.
- **Training is deterministic.** The same CSV and versions reproduce the model file byte for byte.

## Moving to a new scikit-learn version

1. Raise the pin in `pyproject.toml` and run `uv lock`.
2. Run the tests.
   - **If the golden tests pass,** the existing model files still load and predict identically, and you're done.
   - **If they fail,** re-save each bundled model under the new version with the same data and hyperparameters, keeping its id: `scripts/train_primer_scorer.py --training-csv <csv> --id <id> --replace`. Check that the regenerated golden values match the old ones; a mismatch means the new version predicts differently. If they differ, ship the model under a new id instead.

## Training a new model

1. Generate the training CSV from the GAM in qpcr-single-mm, and archive it with the release. Its sha256 is written into the manifest.
2. `uv sync --group train`
3. `uv run --group train python scripts/train_primer_scorer.py --training-csv <csv> --id q5_dct_v2`

   This trains the model, reports its held-out fit, writes `<id>.joblib` and `<id>.json`, and fails unless the bundled model, loaded and encoded as at runtime, reproduces the trained model exactly.
4. To make it the default, update `DEFAULT_MODEL` in `primalscheme3/models/__init__.py`, and add the model to the table above and to `CHANGELOG.md`.
