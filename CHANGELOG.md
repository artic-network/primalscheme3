# Changelog

## 3.5.0 (unreleased, draft)

### Changed
- **primer3-py is no longer used.** Every thermodynamic calculation comes from primalschemers 0.3 (≥0.3, <0.4): thal for hairpins and dimers (identical to primer3 2.3.0), and the same nearest-neighbour model for primer Tm and annealing, so the Python checks and the Rust digest agree.
- **One Tm standard.** The Python checks (downsampling, `downsample-existing-scheme`, reported Tm) now use the digest's: Tm with 0.8 M formamide (new `formamide_conc` config field), failing at or above `primer_tm_max`. Previously they used primer3's Tm with no formamide, about 2 °C higher, and allowed `primer_tm_max + 2`.
- **Primer lengths:** the digest's Tm model changed (primalschemers 0.1 → 0.3), so at about 13% of positions on the test MSA the digest's walk now stops one base shorter at the 5′ end (same 3′ end). Schemes designed with 3.5.0 differ slightly from 3.4.0.
- **The digest gets the run's conditions and hairpin limit.** `mv_conc`, `dv_conc`, `dntp_conc`, `dna_conc`, `formamide_conc` and `primer_hairpin_th_max` are passed to primalschemers instead of its fixed defaults. The hairpin check now runs in Rust (values identical to before); the Python hairpin pass runs only for unchecked digests.
- **Downsampling uses the Rust digest's verdicts.** It digests with `return_all` and takes each sequence's pass/fail from the digest's statuses, instead of re-checking every sequence in Python. Sequences below `min_base_freq`, and positions the digest rejects outright (e.g. Ns without `--ignore-n`), are excluded as before.
- **Downsampled clouds are dimer-checked.** A cloud is rejected if two of its kept primers interact. Previously the downsampling path skipped the dimer check.
- **`repair-mode` reports every reason.** It digests the existing primer sites with `return_all` and logs each sequence's failed checks (e.g. `LowTm, Hairpin`), minor alleles, and how many genomes have no usable sequence and why (e.g. `ContainsN in 2 genomes`).
- **`repair-mode` treats a sequence as already present** when an existing primer shares its 3′ end (one is a 3′-suffix of the other: the same site walked to another length). Previously a prefix match was used.

### Fixed
- `repair-mode` now adds the valid new primers to `primer.bed` (they were reported but never written), and its `config.json` keeps the new MSA's `msa_data` entry.
- `work/file.log` is written by every run, not only the first in a Python process.

### Removed
- The legacy Python digest (`core/digestion.py`, `MSA.digest`); all digestion is in primalschemers.

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
  - `scheme-create` and `panel-create` show `--downsample`, `--downsample-target` and `--primer-scorer-model` in `--help` (they were hidden), and add `--downsample-max-dct` and `--downsample-always-add-prop`.
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
- **Mismatch database log:** `scheme-create` and `panel-create` always logged `Created: …/work/mismatch.db`, though by default the database is held in memory and no file is written. The log now says whether it was built in memory, written to disk (with its path), or not used (`--no-use-matchdb`).
- **Tm reported for an existing bedfile** (`--bedfile` in `scheme-create` / `panel-create`) was about 4 °C too low, because `PrimerPair.calc_tm` passed the DNA concentration (15) as the dNTP concentration (0.8). This only affected that log message: primer selection always used the correct concentrations. Present since 1.3.1 (the June 2024 Typer refactor).

### Removed
- **`core/thermo_misprime.py`**, which could not be imported.
- **`--primer-scorer-model-path`**, replaced by `--primer-scorer-model`. Old `config.json` files that use `primer_scorer_model_path` still load.

### Known issues
These were already present in 3.3.0. They are covered by `expectedFailure` tests and will be fixed when moving to primalschemers 0.2:
- `repair-mode` reports new primers as valid but does not add them to `primer.bed`, and its `config.json` loses the new `msa_data` entry.
- `repair-mode` treats a primer that is a prefix of an existing one as already in the scheme.
- `work/file.log` is only written by the first run in a Python process.
