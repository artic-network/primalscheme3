# Primal-panel

This is a command-line interface tool that generates a primer scheme from a Multiple Sequence Alignment (MSA) file, utalising degenerate primers to handle variation in the genomes.

## Installation

Currently the best way to use is to use poetry to handle dependencies.

```         
git clone https://github.com/ChrisgKent/primal-panel
cd Primal-panel
git submodule add git@github.com:ChrisgKent/primaldimer_py.git
git submodule update --init --recursive
poetry install

poetry run Primal-panel -m msa.fasta -o /path/to/output_dir
```

## Inputs

-   `-m/--msa`: The path to the MSA file(s). If multiple files are to be processed, then pass the paths separated by spaces.
-   `--regionbedfile`: Regions of the MSA that primers should cover.
-   `--primerbedfile`: Add primers to an existing bedfile. Note: The number of pools in bedfile \<= --npools. Primal-panel makes no attempt to validate primers or primer-primer interactions in the bedfile.

## Basic Arguments


-   `-c/--cores`: The number of cores to use in Kmer digestion and thermo checking. Default is `1`.
-   `--ampliconsizemax`: The maximum size of an amplicon. Default is `1000`.
-   `--ampliconsizemin`: The minimum size of an amplicon. Default is `900`.
-   `-o/--output`: The output directory of the primer.bed file. Default is `output`.
-   `--force`: Override the output directory. If set, then any existing output directory will be overwritten.
-   `--minoverlap`: The minimum amount of overlap between primers. Default is `20`.
-   `--npools`: Number of pools to use. Default is `2`.
-   `--minbasefreq`: Primers have have be above this propotion in the input genomes to be included. [`0=<x< 1`]
- `--mapping`: What the primer.bed file uses as a reference. [`first`, `consensus`]
- `--mode`: How amplicons should be selected. [`all`, `region-only`]. `All` picks primers via entropy, `region-only` picks primers for regions spesified in `--regionbedfile`

## Advanced Arguments

-   `--primer_gc_min`: The minimum GC content of a primer. Default is `30`.
-   `--primer_gc_max`: The maximum GC content of a primer. Default is `55`.
-   `--primer_tm_min`: The minimum melting temperature (Tm) of a primer. Default is `59.5`.
-   `--primer_tm_max`: The maximum melting temperature (Tm) of a primer. Default is `62.5`.
-   `--dimerscore`: The threshold for dimer interaction. Default is `-26.0`.
-   `--reducekmers`: Should number of sequences in each Kmer be reduced. Default is `False`.


## Example

```         
poetry run primal-panel -m msa1.fasta msa2.fasta -o /path/to/output_dir --mapping first --minbasefreq 0 --mode all -c 10 --force
```

This command will generate primer schemes for `msa1.fasta` and `msa2.fasta` files and store them in `/path/to/output_dir`. The `--force` option will overwrite any existing output directory.

------------------------------------------------------------------------

This work is licensed under a [Creative Commons Attribution-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-sa/4.0/) 

![](https://i.creativecommons.org/l/by-sa/4.0/88x31.png)
