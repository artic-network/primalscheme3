# Primalscheme3

This is a command-line interface tool that generates a primer scheme from a Multiple Sequence Alignment (MSA) file, utalising degenerate primers to handle variation in the genomes.

## Installation

Currently the best way to use is to use poetry to handle dependencies.

```         
git clone https://github.com/ChrisgKent/primalscheme3
cd primalscheme3
poetry install
poetry build

poetry run primalscheme3 {core-options} {Run Modes} {mode-options}
```

### Core options 

- `-o/--output`: The output directory the scheme in generated into.
- `--force`: Should the output directory be overwritten

- `--primer_gc_min`: Min GC proportion acceptable in primers.
- `--primer_gc_max`: Max GC proportion acceptable in primers.
- `--primer_tm_min`: Min Primer Tm. 
- `--primer_tm_max`: Max Primer Tm.

## Run Modes 
- `scheme-create`: Creates a whole genome overlapping amplicon scheme.
- `scheme-replace`: Finds valid alternative amplicons for an existing amplicon scheme 
- `panel-create`: Creates a amplicon scheme which covers spesified regions of the genomes
- `interaction`: Shows primer-primer interactions located in a primer.bed file

### scheme-create
example: 
```
primalscheme3 -o {output_path} scheme-create --msa {msa_path} --ampliconsize 1000 --mapping first --minbasefreq 0.005 c -10
```

-   `-m/--msa`: The path to the MSA file(s). If multiple files are to be processed, then pass the paths separated by spaces.
-   `-c/--cores`: The number of cores to use in Kmer digestion and thermo checking. Default is `1`.
-   `--ampliconsize`: The size of an amplicon. Either provide one value and amplicons will be within ±10%. Provide two values to set min/max manualy.
-   `--minoverlap`: The minimum amount of overlap between primers. Default is `20`.
-   `--npools`: Number of pools to use. Default is `2`.
-   `--bedfile`: Add primers to an existing bedfile. Note: The number of pools in bedfile \<= --npools. Primal-digest makes no attempt to validate primers or primer-primer interactions in the bedfile.
-   `--mapping`: What should the primers be mapped to. Choice [`consensus`, `first`]. Default is `consensus`.
    -   `consensus`: Uses the MSA consensus as the primer indexing system.
    -   `first`: Uses the first genome in the MSA as the primer indexing system.
-   `--dimerscore`: The threshold for dimer interaction. Default is `-26.0`.
-   `--reducekmers`: Should number of sequences in each Kmer be reduced. Default is `False`. 
-   `--minbasefreq`: The frequency at a SNP/INDEL needs to be at to be included. Default is `0` or all.

### scheme-replace
example: 
```
primalscheme3 -o {output_path} scheme-replace --msa {msa_path} --ampliconsize 1000 --primerbed {path} --config {path}
```

-   `-m/--msa`: The path to the MSA file(s). If multiple files are to be processed, then pass the paths separated by spaces.
-   `--ampliconsize`: The size of an amplicon. Either provide one value and amplicons will be within ±10%. Provide two values to set min/max manualy.
-   `--primerbed`: The primer.bed file of the scheme you want to replace amplicons in.
-   `--ampliconsize`: The size of an amplicon. Either provide one value and amplicons will be within ±10%. Provide two values to set min/max manualy.
-   `--config`: The config.json used to create the original primer scheme

### interactions

```
primalscheme3 -o {anypath} interactions --bedfile {path} --threshold -26
```
-   `--bedfile`: The bedfile to parse for interactions.
-   `--threshold`: Only show dimers below (worse than) this score.

### panel-create

This creates a single pool scheme, in which amplicons cover spesified regions.

#### mode

- ```region-only```: Amplicons only cover region spesified in the --regionbed input
- ```all```: Amplicons cover the regions with the most variance, measured using Shanon Entropy in the MSA.
- `regoion-all`: COMMING SOON. Adds the region amplicons first then keeps adding based on higest entropy

#### --regionbed

This is a bedfile. The 4th colunm is the region name, and the 5th is a score. 
- `chromname`: This needs to match the chromanme used for the msa. Either, first genome ID (--mapping first) or MSA baseame_consensus (--mapping consensus)
- `score`: This is the score an amplicon gets for covering this site. Amplicions with the largest score are added into the scheme first. 

------------------------------------------------------------------------

This work is licensed under a [Creative Commons Attribution-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-sa/4.0/)

![](https://i.creativecommons.org/l/by-sa/4.0/88x31.png)
