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

| H1 | H2 | H3 |
|------------|------|------|
| Test table | col2 | col3 |


`-m  --n-cells`:  Number of cells in pileup file|
` --lambda              required_argument, NULL, 'l'},
{"mu",                  required_argument, NULL, 'u'},
{"p-haploid",           required_argument, NULL, 'h'},
{"p-clonal",            required_argument, NULL, 'c'},
{"pileup-file",         required_argument, NULL, 'p'},
{"vcf-file",            required_argument, NULL, 'o'},
{"n-cells",             required_argument, NULL, 'm'},
{"temp-file",           required_argument, NULL, 'A'},
{"amp-err",             required_argument, NULL, 'B'},
{"p-ado",               required_argument, NULL, 'C'},
{"candidate-threshold", required_argument, NULL, 'D'},
{"posterior-threshold", required_argument, NULL, 'E'},
{"omit-phylo-inference", no_argument     , NULL, 'F'},
{"tree-file",           required_argument, NULL, 'G'}`

Christopher Oldham
University of Connecticut
