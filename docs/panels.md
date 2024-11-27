# Panels

Panels mode was created to enable the targeting of specific regions of the target genome. 

## Modes

There are two main run modes for panels creation.

- `entropy`: Selects regions with the highest entropy

- `region-only`: Selects regions specified in a `.bed` file. See [region-only](###region-only)

### entropy

This mode calculates the shannon entropy of each position. As uses it to determine the most informative region to cover. This approach was designed as a general solution. If specific regions are requires region-only is better.

### region-only

This mode requires a .bed file to define the regions (provided by --region-bedfile)

The region.bed file is a typical .bed file. Specifying the chromosome the region belongs to and its location. 

- `score`: is used to determine the priority of the region when picking amplicons. Ie, a larger score means the region will be added first. 
- `group`: is used to show the regions group. For example, if you have typing scheme for 10 lineages, with each lineages containing 10 snps.  Assigning each SNP to its lineages group allows operations like "Use 2 snps for each group". 
> Lines with either only 5cols or empty (`''`) group col will be treated as normal.



| Col | Meaning     | Type   | Brief description     | Restrictions   |
| --- | ----------- | ------ | --------------------- | -------------- |
| 1   | chrom       | String | Chromosome name       | `[A-Za-z0-9_]` |
| 2   | regionStart | Int    | Region start position | `u64`          |
| 3   | regionEnd   | Int    | Region end position   | `u64`          |
| 4   | regionName  | String | Region name           | `[A-Za-z0-9_]` |
| 5   | score       | Int    | Region score          | `f64`          |
| 6   | group       | String | Region group          | `[A-Za-z0-9_]` |


> The chrom field must match a primary reference in the input msa


