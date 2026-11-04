# Downsampling

> Experimental. Turn it on with `--downsample` on `scheme-create` or `panel-create`.

## Overview

At a variable primer site the genomes don't all share one sequence, so each site has a **primer cloud**: every sequence the genomes have there, with how many genomes have each one. Without downsampling, primalscheme3 orders the whole cloud (every sequence above `--min-base-freq`), and rejects the site if any one of those sequences fails a check (Tm, GC, hairpin, homopolymer, or a dimer within the cloud).

This has two costs on diverse targets:

- **Many primers.** Clouds grow with the number of variants, and every primer is ordered and added to the pool.
- **Lost sites.** One rare variant with a hairpin is enough to lose a site that is otherwise well conserved.

Downsampling relies on a primer still amplifying a template it has a single mismatch with, just a little later. For each cloud it keeps the **smallest useful set of primers** that still amplifies `--downsample-target` of the genomes (99% by default) within `--downsample-max-dct` cycles of a perfect match (−2.0 by default). A primer scorer model predicts the delay (ΔCt) a mismatch causes. See [Primer scorer models](primer-scorer-models.md).

So a rare variant one base away from a common primer doesn't need a primer of its own, and a variant that fails a check no longer rejects the site, as long as another primer covers its genomes.

## The ΔCt model

- **Measured, not assumed.** The model was fitted to qPCR measurements of 270 single-mismatch primers (SARS-CoV-2 targets, Q5 polymerase), each compared with its perfectly matched primer on the same plate.
- **What it uses.** The prediction comes from the stability of the primer beyond the mismatch, how far the mismatch is from the 3′ end, the substitution, and its 3′ neighbour. For Q5, the mismatch's position and its downstream context matter more than how much it weakens binding.
- **Yes or no at −2.** A ΔCt of −2 is about a 4-fold drop in yield. The prediction is only used to decide whether a primer covers a template or not. It isn't used to rank primers.
- **Accuracy.** Predictions are typically within about 1.4 cycles of the measured value. As a yes/no call at −2 the model is reliable (ROC AUC 0.93).
- **Cautious by design.** A mismatch outside the range the model was trained on gets no prediction, so it doesn't count as covered. Predictions from several primers aren't added together. Both rules make it keep more primers rather than fewer.

See [Primer scorer models](primer-scorer-models.md) for how the model is built and bundled.

## How a cloud is downsampled

### 1. The cloud

The digest returns every sequence at the site, with its genome count and which checks it fails. Sequences below `--min-base-freq` are dropped, as in a normal digest. Sequences that fail a check, including one that dimerises with itself, **can't be chosen**, but their genomes still need covering.

### 2. Who covers whom

Two sequences are **neighbours** if they differ at exactly one position, compared from the 3′ end. For each pair of neighbours, the model predicts the ΔCt of each primer on the other's template, both ways round.

A primer **covers** a genome if the genome's sequence:

- is the primer's own sequence (ΔCt 0), or
- is a neighbour, and the predicted ΔCt is at least `--downsample-max-dct`.

A genome counts as covered by the best kept primer for it. Genomes two or more mismatches from every kept primer are not covered.

For example, a cloud at one site across 100 genomes. The ΔCt values are `q5_dct_v1`'s predictions for the first primer on each template. The pass/fail statuses are made up for the example.

```
  sequence              genomes  checks          ΔCt of primer 1   covered by primer 1
1 ACGTTGCAAGTCCGATGCA   80       passes           0     (itself)    yes
2 ACGTTGCAATTCCGATGCA   10       passes          -0.44             yes
3 ACGTTGCAAGTCCGCTGCA    5       fails Hairpin   -0.69             yes
4 ACGTTGCAAGTCCGATGAA    4       passes          -4.86             no
5 ACGTAGCAAGTCCGTTGCA    1       passes          (two mismatches)  no
```

- Primer 1 is in 80% of genomes, so it is always kept, and covers 95 genomes: its own and those of sequences 2 and 3. Sequence 3 fails a check and can't be ordered, but its genomes are still amplified.
- 95% is below the 99% target, so the greedy step adds sequence 4. Its mismatch is next to the 3′ end and predicted to cost almost five cycles, so it needs a primer of its own. That brings coverage to 99%.
- The cloud is kept with two primers instead of four orderable ones. The site would otherwise have been rejected outright because of sequence 3. The one genome with sequence 5 isn't covered, which the target allows.

How much a mismatch costs depends on its position and the substitution. For single mismatches across this primer, the model predicts from +0.5 to −4.9, so whether a variant needs its own primer is decided per variant, not by a rule of thumb.

### 3. Keep the common sequences

Every passing sequence found in at least `--downsample-always-add-prop` of the genomes (25% by default) is kept, whether or not it is needed for coverage.

### 4. Add primers greedily

Until `--downsample-target` of the genomes are covered, it adds the passing primer that covers the most genomes not yet covered. Ties go to the more common sequence. A primer that dimerises with one already kept is skipped in favour of the next best. It stops early if no primer adds coverage.

### 5. The verdict

The site is kept, with just the chosen primers, if the target is reached. It is rejected if the target can't be reached, or if two of the always-kept primers dimerise with each other. The log records why, for example `BelowTarget`, or the checks that stopped the uncovered genomes' sequences from being chosen (such as `Hairpin`).

## Not an optimal search

Downsampling is a fast per-site heuristic. It doesn't find the best possible primer set or model every way a primer could amplify. In particular:

- **Greedy, not minimal.** Each step takes the primer that helps most right then, which usually gives a small set but isn't guaranteed to give the smallest.
- **Single mismatches only.** Only neighbours one mismatch apart are scored. A primer is never credited with amplifying a template two or more mismatches away, even if in practice it would.
- **One model, one polymerase.** ΔCt comes from `q5_dct_v1`, which was trained on Q5 qPCR data for a single mismatch. Other polymerases, conditions and mismatch combinations aren't modelled.
- **Best primer only.** A genome counts as covered by its single best primer. Contributions from several primers aren't added up.
- **Genome counts, not prevalence.** Coverage is weighted by how many genomes in the MSA have each sequence, so it reflects how the MSA was sampled.
- **One site at a time.** Downsampling only looks at its own cloud. Dimers with other primers in the pool, and mispriming, are checked later, when the pair is added to a pool.

## Using it

```bash
primalscheme3 scheme-create --msa target.fasta --output out --downsample
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--downsample` | off | Turn downsampling on |
| `--downsample-target` | 0.99 | Proportion of genomes each cloud must cover |
| `--downsample-max-dct` | −2.0 | Lowest predicted ΔCt at which a primer still covers a genome |
| `--downsample-always-add-prop` | 0.25 | Always keep passing primers found in at least this proportion of genomes |
| `--primer-scorer-model` | `q5_dct_v1` | Bundled model id, or the path to a model's `.json` manifest |

Downsampling runs in parallel with `--ncores`. The model used and its sha256 are recorded in `config.json` (`primer_scorer_model`, `primer_scorer_model_sha256`).

With the `downsample_plot` config flag set, each cloud's neighbour graph is drawn as an SVG, and a `downsample_summary.csv` of each cloud's size before and after downsampling and its coverage is written, both under `work/downsample_plots`. This is slow, so it's meant for checking individual sites.
