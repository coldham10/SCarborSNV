# SCarborSNV

**SCarborSNV** is an efficient algorithm for phylogeny aware genotyping of single cell genomes.
**SCarborSNV** uses dynamic programming to compute the expected evolutionary distance between pairs of cells, then missing information is inferred from an approximate phylogeny built by neighbor joining.

## Installation

Clone **SCarborSNV** from github:
`git clone https://github.com/coldham10/SCarborSNV.git`

To install, in the newly created directory:
`make`

## Running **ScarborSNV**

**SCarborSNV** takes as input an mpileup file or stream. You must specify the number of cells in the pieup with the `-m` option.

From a pre-existing mpileup file with 10 cells:
`./SCarborSNV -m 10 -p tencells.mpileup`

Or, for example, to use a compressed mpileup file:
`zcat twentycells.mpileup.gz | ./SCarborSNV -m 20`

Command line options are:

| Short option | Long option | Description |
|------------|------|------|
| `-m` | `--n-cells` | (Required) Number of cells in the mpileup |
| `-p` | `--pileup-file` | Name of mpileup file to read |
| `-o` | `--vcf-file` | Name of output file. Defaults to SCarborSNV\_out.vcf |
|      | `--lambda` | Somatic mutation rate |
|      | `--mu` | Germline mutation rate |
|      | `--p-haploid` | Prior probability that any given locus has experienced LOH |
|      | `--p-clonal` | Prior probability that any SNV is public |
|      | `--temp-file` | Name of necessary temp file. By default is written to /tmp. **Warning**: If running multiple instances of **SCarborSNV** concurrently you **must** specify different temp file locations |
|      | `--amp-err` | Probability of a base substitution due to amplification error |
|      | `--p-ado` | Prior probability that a base has experienced allelic dropout |
|  | `--candidate-threshold` | Threshold on posterior probability of alternate allele count being 0 at a locus. Loci with P(aac = 0) < threshold will be considered as candidate mutants. |

## Licence
The **SCarborSNV** software is freely available under the GPL v.3. 

Copyright 2019 Christopher Oldham




Christopher Oldham

University of Connecticut
