# Primalscheme3

[![CI](https://github.com/ChrisgKent/primalscheme3/actions/workflows/pytest.yml/badge.svg)](https://github.com/ChrisgKent/primalscheme3/actions/workflows/pytest.yml) [![Generic badge](https://img.shields.io/badge/DOI-10.1101/2024.12.20.629611-blue.svg)](https://doi.org/10.1101/2024.12.20.629611)

This is a command-line interface tool that generates a primer scheme from a Multiple Sequence Alignment (MSA) file, utilising degenerate primers to handle variation in the genomes.

## Installation

Via pip!
```
pip install primalscheme3
```

From source 
```         
git clone https://github.com/artic-network/primalscheme3
cd primalscheme3
uv sync
uv run primalscheme3
```

# `PrimalScheme3`

**Usage**:

```console
$ primalscheme3 [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--version`
* `--install-completion`: Install completion for the current shell.
* `--show-completion`: Show completion for the current shell, to copy it or customize the installation.
* `--help`: Show this message and exit.

**Commands**:

* `scheme-create`: Creates a tiling overlap scheme for each...
* `replace-primerpair`: Replaces a primerpair in a bedfile
* `panel-create`: Creates a primer panel
* `interactions`: Shows all the primer-primer interactions...
* `repair-mode`: Repairs a primer scheme via adding more...
* `visualise-primer-mismatches`: Visualise mismatches between primers and...
* `visualise-primer-dct`: Visualise primer mismatch penalty (DCT)...
* `visualise-bedfile`: Visualise the bedfile

## `primalscheme3 scheme-create`

Creates a tiling overlap scheme for each MSA file

**Usage**:

```console
$ primalscheme3 scheme-create [OPTIONS]
```

**Options**:

* `--msa PATH`: The MSA to design against. To use multiple MSAs, use multiple --msa flags. (--msa 1.fasta --msa 2.fasta)  [required]
* `--output PATH`: The output directory  [required]
* `--amplicon-size INTEGER RANGE`: The size of an amplicon. Min / max size are ± 10 percent [100&lt;=x&lt;=2000]  [default: 400; 100&lt;=x&lt;=2000]
* `--bedfile PATH`: An existing bedfile to add primers to
* `--min-overlap INTEGER RANGE`: min amount of overlap between primers  [default: 10; x&gt;=0]
* `--n-pools INTEGER RANGE`: Number of pools to use  [default: 2; x&gt;=1]
* `--dimer-score FLOAT`: Threshold for dimer interaction  [default: -26.0]
* `--min-base-freq FLOAT RANGE`: Min freq to be included,[0&lt;=x&lt;=1]  [default: 0.0; 0.0&lt;=x&lt;=1.0]
* `--mapping [first|consensus]`: How should the primers in the bedfile be mapped  [default: first]
* `--circular / --no-circular`: Should a circular amplicon be added  [default: no-circular]
* `--backtrack / --no-backtrack`: Should the algorithm backtrack  [default: no-backtrack]
* `--ignore-n / --no-ignore-n`: Should N in the input genomes be ignored  [default: no-ignore-n]
* `--force / --no-force`: Override the output directory  [default: no-force]
* `--input-bedfile PATH`: Path to a primer.bedfile containing the pre-calculated primers
* `--high-gc / --no-high-gc`: Use high GC primers  [default: no-high-gc]
* `--offline-plots / --online-plots`: Offline plots includes 3Mb of dependencies, so they can be viewed offline  [default: offline-plots]
* `--use-matchdb / --no-use-matchdb`: Create and use a mispriming database  [default: use-matchdb]
* `--ncores INTEGER RANGE`: Number of CPU cores to use during digestion  [default: 1; x&gt;=1]
* `--use-annealing / --use-tm`: Using annealing proportion rather than Tm to calculate primers  [default: use-tm]
* `--help`: Show this message and exit.

## `primalscheme3 replace-primerpair`

Replaces a primerpair in a bedfile

**Usage**:

```console
$ primalscheme3 replace-primerpair [OPTIONS] PRIMERNAME PRIMERBED MSA
```

**Arguments**:

* `PRIMERNAME`: The name of the primer to replace  [required]
* `PRIMERBED`: The bedfile containing the primer to replace  [required]
* `MSA`: The msa used to create the original primer scheme  [required]

**Options**:

* `--amplicon-size-max INTEGER RANGE`: The max size of an amplicon [100&lt;=x&lt;=2000]  [100&lt;=x&lt;=2000; required]
* `--amplicon-size-min INTEGER RANGE`: The min size of an amplicon. [100&lt;=x&lt;=2000]  [100&lt;=x&lt;=2000; required]
* `--config-path PATH`: The config.json used to create the original primer scheme  [required]
* `--output PATH`: The output directory  [required]
* `--mode [addbest|listall]`: Select what run mode  [default: listall]
* `--force / --no-force`: Override the output directory  [default: no-force]
* `--mask-old-sites / --no-mask-old-sites`: If True prevents replacement primers from spanning old primer regions.  [default: mask-old-sites]
* `--help`: Show this message and exit.

## `primalscheme3 panel-create`

Creates a primer panel

**Usage**:

```console
$ primalscheme3 panel-create [OPTIONS]
```

**Options**:

* `--msa PATH`: Paths to the MSA files  [required]
* `--output PATH`: The output directory  [required]
* `--region-bedfile FILE`: Path to the bedfile containing the wanted regions
* `--input-bedfile FILE`: Path to a primer.bedfile containing the pre-calculated primers
* `--mode [entropy|region-only|equal]`: Select what run mode  [default: region-only]
* `--amplicon-size INTEGER`: The size of an amplicon  [default: 400]
* `--n-pools INTEGER RANGE`: Number of pools to use  [default: 2; x&gt;=1]
* `--dimer-score FLOAT`: Threshold for dimer interaction  [default: -26.0]
* `--min-base-freq FLOAT RANGE`: Min freq to be included,[0&lt;=x&lt;=1]  [default: 0.0; 0.0&lt;=x&lt;=1.0]
* `--mapping [first|consensus]`: How should the primers in the bedfile be mapped  [default: first]
* `--max-amplicons INTEGER RANGE`: Max number of amplicons to create  [x&gt;=1]
* `--max-amplicons-msa INTEGER RANGE`: Max number of amplicons for each MSA  [x&gt;=1]
* `--max-amplicons-region-group INTEGER RANGE`: Max number of amplicons for each region  [x&gt;=1]
* `--force / --no-force`: Override the output directory  [default: no-force]
* `--high-gc / --no-high-gc`: Use high GC primers  [default: no-high-gc]
* `--offline-plots / --online-plots`: Offline plots includes 3Mb of dependencies, so they can be viewed offline  [default: offline-plots]
* `--use-matchdb / --no-use-matchdb`: Create and use a mispriming database  [default: use-matchdb]
* `--ncores INTEGER RANGE`: Number of CPU cores to use during digestion  [default: 1; x&gt;=1]
* `--use-annealing / --use-tm`: Using annealing proportion rather than Tm to calculate primers  [default: use-tm]
* `--help`: Show this message and exit.

## `primalscheme3 interactions`

Shows all the primer-primer interactions within a bedfile

**Usage**:

```console
$ primalscheme3 interactions [OPTIONS] BEDFILE
```

**Arguments**:

* `BEDFILE`: Path to the bedfile  [required]

**Options**:

* `--threshold FLOAT`: Only show interactions more severe (Lower score) than this value  [default: -26.0]
* `--help`: Show this message and exit.

## `primalscheme3 repair-mode`

Repairs a primer scheme via adding more primers to account for new mutations

**Usage**:

```console
$ primalscheme3 repair-mode [OPTIONS]
```

**Options**:

* `--bedfile PATH`: Path to the bedfile  [required]
* `--msa PATH`: An MSA, with the reference.fasta, aligned to any new genomes with mutations  [required]
* `--config PATH`: Path to the config.json  [required]
* `--output PATH`: The output directory  [required]
* `--force / --no-force`: Override the output directory  [default: no-force]
* `--help`: Show this message and exit.

## `primalscheme3 visualise-primer-mismatches`

Visualise mismatches between primers and the input genomes

**Usage**:

```console
$ primalscheme3 visualise-primer-mismatches [OPTIONS] MSA BEDFILE
```

**Arguments**:

* `MSA`: The MSA used to design the scheme  [required]
* `BEDFILE`: The bedfile containing the primers  [required]

**Options**:

* `--output FILE`: Output location of the plot  [default: primer.html]
* `--include-seqs / --no-include-seqs`: Reduces plot filesize, by excluding primer sequences  [default: include-seqs]
* `--offline-plots / --no-offline-plots`: Includes 3Mb of dependencies into the plots, so they can be viewed offline  [default: offline-plots]
* `--help`: Show this message and exit.

## `primalscheme3 visualise-primer-dct`

Visualise primer mismatch penalty (DCT) between primers and the input genomes

**Usage**:

```console
$ primalscheme3 visualise-primer-dct [OPTIONS] MSA BEDFILE
```

**Arguments**:

* `MSA`: The MSA used to design the scheme  [required]
* `BEDFILE`: The bedfile containing the primers  [required]

**Options**:

* `--model TEXT`: Primer scorer model: a bundled model id, or the path to a model&#x27;s .json manifest  [default: q5_dct_v1]
* `--output FILE`: Output location of the plot  [default: primer_dct.html]
* `--include-seqs / --no-include-seqs`: Reduces plot filesize, by excluding primer sequences  [default: include-seqs]
* `--offline-plots / --no-offline-plots`: Includes 3Mb of dependencies into the plots, so they can be viewed offline  [default: offline-plots]
* `--help`: Show this message and exit.

## `primalscheme3 visualise-bedfile`

Visualise the bedfile

**Usage**:

```console
$ primalscheme3 visualise-bedfile [OPTIONS] BEDFILE REF_PATH
```

**Arguments**:

* `BEDFILE`: The bedfile containing the primers  [required]
* `REF_PATH`: The reference.fasta file containing the --ref-id genome  [required]

**Options**:

* `--ref-id TEXT`: The reference genome ID  [required]
* `--output FILE`: Output location of the plot  [default: bedfile.html]
* `--help`: Show this message and exit.
