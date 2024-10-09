# Primalscheme3

[![CI](https://github.com/ChrisgKent/primalscheme3/actions/workflows/pytest.yml/badge.svg)](https://github.com/ChrisgKent/primalscheme3/actions/workflows/pytest.yml)

This is a command-line interface tool that generates a primer scheme from a Multiple Sequence Alignment (MSA) file, utalising degenerate primers to handle variation in the genomes.

## Installation

Currently the best way to use is to use poetry to handle dependencies.

```         
git clone https://github.com/ChrisgKent/primalscheme3
cd primalscheme3
poetry install
poetry build

```

# `PrimalScheme3`

**Usage**:

```console
$ PrimalScheme3 [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--version`
* `--install-completion`: Install completion for the current shell.
* `--show-completion`: Show completion for the current shell, to copy it or customize the installation.
* `--help`: Show this message and exit.

**Commands**:

* `interactions`: Shows all the primer-primer interactions...
* `panel-create`: Creates a primerpanel
* `remap-mode`: Remaps a primer scheme to a new reference...
* `repair-mode`: Repairs a primer scheme via adding more...
* `scheme-create`: Creates a tiling overlap scheme for each...
* `scheme-replace`: Replaces a primerpair in a bedfile

## `PrimalScheme3 interactions`

Shows all the primer-primer interactions within a bedfile

**Usage**:

```console
$ PrimalScheme3 interactions [OPTIONS] BEDFILE
```

**Arguments**:

* `BEDFILE`: Path to the bedfile  [required]

**Options**:

* `--threshold FLOAT`: Only show interactions more severe (Lower score) than this value  [default: -26.0]
* `--help`: Show this message and exit.

## `PrimalScheme3 panel-create`

Creates a primerpanel

**Usage**:

```console
$ PrimalScheme3 panel-create [OPTIONS]
```

**Options**:

* `--msa PATH`: Paths to the MSA files  [required]
* `--output PATH`: The output directory  [required]
* `--regionbedfile PATH`: Path to the bedfile containing the wanted regions
* `--inputbedfile PATH`: Path to a primer.bedfile containing the precalculated primers
* `--mode [all|region-only|region-all]`: Select what mode for selecting regions in --regionbedfile  [default: region-only]
* `--amplicon-size INTEGER`: The size of an amplicon  [default: 400]
* `--n-pools INTEGER RANGE`: Number of pools to use  [default: 2; x>=1]
* `--dimer-score FLOAT`: Threshold for dimer interaction  [default: -26.0]
* `--min-base-freq FLOAT RANGE`: Min freq to be included,[0<=x<=1]  [default: 0.0; 0.0<=x<=1.0]
* `--mapping [first|consensus]`: How should the primers in the bedfile be mapped  [default: first]
* `--maxamplicons INTEGER RANGE`: Max number of amplicons to create  [x>=1]
* `--force / --no-force`: Override the output directory  [default: no-force]
* `--high-gc / --no-high-gc`: Use high GC primers  [default: no-high-gc]
* `--help`: Show this message and exit.

## `PrimalScheme3 remap-mode`

Remaps a primer scheme to a new reference genome

**Usage**:

```console
$ PrimalScheme3 remap-mode [OPTIONS]
```

**Options**:

* `--bedfile PATH`: Path to the bedfile  [required]
* `--id-to-remap-to TEXT`: The ID of the reference genome to remap to  [required]
* `--msa PATH`: Path to the MSA file  [required]
* `--output PATH`: The output directory  [required]
* `--help`: Show this message and exit.

## `PrimalScheme3 repair-mode`

Repairs a primer scheme via adding more primers to account for new mutations

**Usage**:

```console
$ PrimalScheme3 repair-mode [OPTIONS]
```

**Options**:

* `--bedfile PATH`: Path to the bedfile  [required]
* `--msa PATH`: An MSA, with the reference.fasta, aligned to any new genomes with mutations  [required]
* `--config PATH`: Path to the config.json  [required]
* `--output PATH`: The output directory  [required]
* `--force / --no-force`: Override the output directory  [default: no-force]
* `--help`: Show this message and exit.

## `PrimalScheme3 scheme-create`

Creates a tiling overlap scheme for each MSA file

**Usage**:

```console
$ PrimalScheme3 scheme-create [OPTIONS]
```

**Options**:

* `--msa PATH`: The name of the scheme  [required]
* `--output PATH`: The output directory  [required]
* `--amplicon-size INTEGER`: The size of an amplicon. Use single value for ± 10 percent [100<=x<=2000]  [default: 400]
* `--bedfile PATH`: An existing bedfile to add primers to
* `--min-overlap INTEGER RANGE`: min amount of overlap between primers  [default: 10; x>=0]
* `--n-pools INTEGER RANGE`: Number of pools to use  [default: 2; x>=1]
* `--dimer-score FLOAT`: Threshold for dimer interaction  [default: -26.0]
* `--min-base-freq FLOAT RANGE`: Min freq to be included,[0<=x<=1]  [default: 0.0; 0.0<=x<=1.0]
* `--mapping [first|consensus]`: How should the primers in the bedfile be mapped  [default: first]
* `--circular / --no-circular`: Should a circular amplicon be added  [default: no-circular]
* `--backtrack / --no-backtrack`: Should the algorithm backtrack  [default: no-backtrack]
* `--ignore-n / --no-ignore-n`: Should N in the input genomes be ignored  [default: no-ignore-n]
* `--force / --no-force`: Override the output directory  [default: no-force]
* `--input-bedfile PATH`: Path to a primer.bedfile containing the precalculated primers
* `--high-gc / --no-high-gc`: Use high GC primers  [default: no-high-gc]
* `--help`: Show this message and exit.

## `PrimalScheme3 scheme-replace`

Replaces a primerpair in a bedfile

**Usage**:

```console
$ PrimalScheme3 scheme-replace [OPTIONS] PRIMERNAME PRIMERBED MSA
```

**Arguments**:

* `PRIMERNAME`: The name of the primer to replace  [required]
* `PRIMERBED`: The bedfile containing the primer to replace  [required]
* `MSA`: The msa used to create the original primer scheme  [required]

**Options**:

* `--ampliconsize INTEGER`: The size of an amplicon. Use single value for ± 10 percent [100<=x<=2000]  [required]
* `--config PATH`: The config.json used to create the original primer scheme  [required]
* `--help`: Show this message and exit.

------------------------------------------------------------------------

This work is licensed under GNU General Public License v3.
