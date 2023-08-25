# Primal-digest

This is a command-line interface tool that generates a primer scheme from a Multiple Sequence Alignment (MSA) file, utalising degenerate primers to handle variation in the genomes.

## Installation

Currently the best way to use is to use poetry to handle dependencies.

```         
git clone https://github.com/ChrisgKent/primal-digest
cd primal-digest
poetry install

poetry run primal-digest -m msa.fasta -o /path/to/output_dir
```

## Basic Arguments

-   `-m/--msa`: The path to the MSA file(s). If multiple files are to be processed, then pass the paths separated by spaces.
-   `-c/--cores`: The number of cores to use in Kmer digestion and thermo checking. Default is `1`.
-   `--ampliconsizemax`: The maximum size of an amplicon. Default is `1000`.
-   `--ampliconsizemin`: The minimum size of an amplicon. Default is `900`.
-   `-p/--prefix`: The prefix used in the bedfile name. Default is `output`.
-   `-o/--output`: The output directory of the primer.bed file. Default is `output`.
-   `--refnames`: The name of the reference used in the outputted bed file. If processing multiple MSA files, then pass the reference names separated by spaces. Default is `MSA`.
-   `--force`: Override the output directory. If set, then any existing output directory will be overwritten.
-   `--minoverlap`: The minimum amount of overlap between primers. Default is `20`.
-   `--npools`: Number of pools to use. Default is `2`.
-   `--bedfile`: Add primers to an existing bedfile. Note: The number of pools in bedfile \<= --npools. Primal-digest makes no attempt to validate primers or primer-primer interactions in the bedfile.

## Advanced Arguments

-   `--primer_gc_min`: The minimum GC content of a primer. Default is `30`.
-   `--primer_gc_max`: The maximum GC content of a primer. Default is `55`.
-   `--primer_tm_min`: The minimum melting temperature (Tm) of a primer. Default is `59.5`.
-   `--primer_tm_max`: The maximum melting temperature (Tm) of a primer. Default is `62.5`.
-   `--dimerscore`: The threshold for dimer interaction. Default is `-26.0`.
-   `--reducekmers`: Should number of sequences in each Kmer be reduced. Default is `False`.
-   `--minbasefreq`: The frequency at a SNP/INDEL needs to be at to be included. Default is `0` or all.




## Example

```         
poetry run primal-digest -m msa1.fasta msa2.fasta -o /path/to/output_dir --refnames MSA1 MSA2 --force
```

This command will generate primer schemes for `msa1.fasta` and `msa2.fasta` files and store them in `/path/to/output_dir`. The `--force` option will overwrite any existing output directory.

------------------------------------------------------------------------

This work is licensed under a [Creative Commons Attribution-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-sa/4.0/) 

![](https://i.creativecommons.org/l/by-sa/4.0/88x31.png)
