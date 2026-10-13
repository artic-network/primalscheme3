# Changelog

## 3.4.0 (unreleased, draft)

### Added
- **Experimental ΔCt downsampling** (`scheme-create` / `panel-create` `--downsample`). This reduces the number of primers in each degenerate primer cloud. It keeps the primers needed so that `--downsample-target` of the genomes are amplified within `downsample_max_dct` (default −2.0) of a perfect match, as predicted by a primer scorer model. Primers that fail thermodynamic checks are never kept.
- **Bundled primer scorer model `q5_dct_v1`**. It predicts the ΔCt of a primer with one mismatch to its template (Q5 polymerase; a scikit-learn `HistGradientBoostingRegressor` surrogate of a GAM fitted to qPCR data; held-out fit R² 0.99957).
  - It ships inside the package, so no model file or download is needed.
  - Choose a model with `--primer-scorer-model <id or manifest path>` or `PRIMALSCHEME3_PRIMER_SCORER_MODEL`.
  - Features are always computed at the conditions the model was trained at, whatever the run's settings.
  - The model id and sha256 used are recorded in `config.json` (`primer_scorer_model`, `primer_scorer_model_sha256`).
  - Its manifest records the feature encoding, the training data's sha256, the fit, and golden predictions that CI checks.
  - See `docs/primer-scorer-models.md` for how the model is made, the versioning rules, and how to train a new model.
- **`visualise-primer-dct`**: a heatmap of the predicted mismatch penalty (ΔCt) between each primer and each genome in an MSA.
- **Downsampling options and output**:
  - The downsampling step runs in parallel with `--ncores`.
  - The optional `downsample_plot` config flag writes a per-cloud graph (SVG) and `downsample_summary.csv` to `work/downsample_plots`.
- **`ps3` command**, a short alias for `primalscheme3`.
- **Docker image** built with uv (multi-stage, runs as a non-root user).
- **`scripts/train_primer_scorer.py`**, which trains a primer scorer model, reports its fit, and bundles it with its manifest. Retraining from the same data and versions is byte-for-byte reproducible (`uv sync --group train`).

### Changed
- **`downsample-existing-scheme`** (hidden) now uses the ΔCt model, the same as `scheme-create --downsample`, in place of annealing. It keeps the bed header and probe lines, and handles amplicons with only LEFT or only RIGHT primers. It also takes `--primer-scorer-model` and `--downsample-max-dct`.
- **Chrom names:** every character outside primalbedtools' `CHROM_REGEX` is now replaced with `_`, not just `/` and `-`.
- **`scheme-replace`:** replacement primers are sorted by coverage, and it now writes `config.json` and the amplicon bed files.
- **Remapping:** `BedPrimerPair`s are remapped on genome-mapped MSAs.
- **Dependencies:**
  - primalschemers ≥ 0.1.13 is required.
  - New runtime dependencies: `scikit-learn` (pinned to 1.8.x, the version the model was saved with), `joblib`, `threadpoolctl` and `scipy`. scikit-learn's wheels bundle their own OpenMP, so nothing extra is needed on macOS.
  - `pandas` and `pyarrow` are needed only for training (`train` dependency group).
- The CLI reference (`docs/cli.md` and the README) has been regenerated.

### Fixed
- **Downsampling:**
  - primers of different lengths in one cloud were not compared;
  - mismatch scoring failed when primers differed in length (`.` padding).
- **Repair mode:** primer pools are sized by the highest pool index, so sparse pool numbering no longer raises an IndexError.
- **`check_for_end_on_gap`:** guarded against index 0.
- **Tm reported for an existing bedfile** (`--bedfile` in `scheme-create` / `panel-create`) was about 4 °C too low, because `PrimerPair.calc_tm` passed the DNA concentration (15) as the dNTP concentration (0.8). This only affected that log message: primer selection always used the correct concentrations. Present since 1.3.1 (the June 2024 Typer refactor).

### Removed
- **`core/thermo_misprime.py`**, which could not be imported.
- **`--primer-scorer-model-path`**, replaced by `--primer-scorer-model`. Old `config.json` files that use `primer_scorer_model_path` still load.

### Known issues
These were already present in 3.3.0. They are covered by `expectedFailure` tests and will be fixed when moving to primalschemers 0.2:
- `repair-mode` reports new primers as valid but does not add them to `primer.bed`, and its `config.json` loses the new `msa_data` entry.
- `repair-mode` treats a primer that is a prefix of an existing one as already in the scheme.
- `work/file.log` is only written by the first run in a Python process.
