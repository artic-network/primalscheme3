# Mispriming

## Overview

A primer can bind places other than the site it was designed for: elsewhere in the target genome, in another MSA in the same run, or in a background genome such as the host. A single stray binding site is usually harmless, because one primer on its own only copies a template linearly. The problem is a **product**: a forward-facing site and a reverse-facing site on the same sequence, close enough together to be amplified. Off-target products compete with the real amplicons and give off-target reads.

The mispriming check stops these products. When primalscheme3 tries to add a primer pair to a pool, it rejects the pair if doing so would create an off-target product shorter than `mismatch_product_size` (by default `amplicon_size_max`), either:

- **with the pool**: one of the pair's binding sites and a site of a primer already in the pool, or
- **alone**: two of the pair's own binding sites, other than its intended amplicon.

It runs by default in `scheme-create` and `panel-create`. `--no-use-matchdb` turns it off.

## Finding binding sites

### Only the 3′ end

Polymerase extends a primer from its 3′ end. A mismatch near the 3′ end stops extension, while the 5′ end tolerates mismatches. So the check uses only each primer's **3′ terminal kmer**: its last `k` bases, where `k` is `primer_size_min` (19 by default, 17 with `--high-gc`).

A site that matches those `k` bases is a binding site, whatever the rest of the primer does there.

### One mismatch

A site counts if it matches the 3′ kmer **exactly or with one substitution** (`mismatch_fuzzy`, on by default). In other words, it allows a Hamming distance of at most 1 over the kmer.

Every kmer is expanded into itself plus its `3k` single-substitution variants (57 for a 19-mer). All of them go into one index, and a reference window is a hit if it is in that index. A 19-base window that differs from the primer's 3′ end at one position is found. A window that differs at two positions, or has an insertion or deletion, is not.

```
primer 3′ kmer      ACGTTGCAAGTCCGATGCA
reference window    ACGTTGCAAGTACGATGCA   one mismatch    → binding site
reference window    ACGTTGCAAGTACGATGGA   two mismatches  → not a binding site
```

### Every reference, both strands

The kmers are looked for on both strands of:

- every genome in every MSA in the run, with alignment gaps removed, and
- every record in each `--background` FASTA (plain or gzipped, any size).

A degenerate primer is screened as each sequence it expands to. Ambiguity codes in a reference window are expanded too. Windows that would expand to more than 16 sequences, or that contain an `N`, are skipped.

### Built for large genomes

The primers are indexed, not the genomes. All candidate primers are digested first, then their 3′ kmers and variants are put into a single index, and each reference is streamed past it once. Memory depends on the number of primers rather than the size of the references, so a host genome can be screened as a background without building a database of its kmers. Hits are cached by kmer, so each one is only looked up once per run.

## Finding products

Each binding site is recorded as a hit:

| Field | Meaning |
| --- | --- |
| reference | Which MSA or background file |
| contig | Which record in a background file (always 0 for an MSA) |
| position | Where the primer's 3′ end base sits |
| strand | `+` if the primer points right, `-` if it points left |

Two hits form a product when they are on the same contig, the `+` hit is upstream of the `-` hit, and they are less than `mismatch_product_size` apart:

```
        + hit                                  - hit
  5′ ───►───────────────────────────────────────◄─── 3′
        |<────── less than mismatch_product_size ─────>|   → product
```

Each pool keeps every binding site of every primer in it, including their intended sites. When a pair is tried in a pool, the check looks for a product:

1. **With the pool**: between any of the new pair's hits and the pool's sites.
2. **Alone**: between the pair's own hits. The pair's two intended sites (its real amplicon) are allowed to pair with each other, and any other combination counts.

The first product found rejects the pair for that pool, and the reference it is on is recorded.

## Not an exhaustive search

The check is a fast filter for near-identical 3′ ends against defined background. It doesn't model every place a primer could bind, and there's no alignment or thermodynamic scoring of off-target sites. In particular:

- **Two or more mismatches** in the 3′ kmer, or an **insertion or deletion between the primer and the site**, are not found. Alignment gaps are different: they are removed before scanning, so a site spanning a gap in the MSA is found.
- **The 5′ end is ignored.** A site matching the last `k` bases is a hit even if the rest of the primer doesn't match there.
- **Mismatch position is not weighted.** A single mismatch counts wherever it is in the kmer, even at the 3′ base where extension is unlikely. This errs on the side of rejecting a pair.
- **Products longer than `mismatch_product_size`** are ignored.
- **Primer-dimers** are a separate check (`--dimer-score`).

## Using it

```bash
primalscheme3 scheme-create --msa target.fasta --output out \
    --background host.fasta.gz --background other.fasta
```

- `--background PATH` adds a genome to screen against. It can be given more than once. Backgrounds are ignored with `--no-use-matchdb`.
- `mismatch_fuzzy`, `mismatch_kmersize` (set from `primer_size_min`) and `mismatch_product_size` are recorded in `config.json`.

### Output

The run log shows how many primer pair checks passed or failed at each stage, and which references the rejected products were on:

- **Primer pair checks**: how many pairs were rejected for overlap, interaction or mispriming, or were OK. A pair is counted each time it is tried in a pool, and stops at its first failing check.
- **Mispriming products by reference**: for each MSA or background, how many products were formed by a pair **alone** or **with the pool**.

The same counts go into `config.json` under `pair_checks`. Background paths go under `background_genomes`, and each file's md5 under `background_data`.
