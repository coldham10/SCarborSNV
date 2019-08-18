/**
 * SCarborSNV: Efficient Phylogeny-aware Single Nucleotide Variant Detection for Single Cells
 *
 * Copyright (C) 2019 Christopher Oldham
 *
 * This file is part of SCarborSNV.
 *
 * SCarborSNV is free software: you can redistribute it and/or modify
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SCarborSNV is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SCarborSNV.  If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef SEQUENCE_UTILS_H_
#define SEQUENCE_UTILS_H_

#include <math.h>
#include "scarborsnv.h"
#include "pileup_reader.h"


/*Multiply a phred by this to get the log probability of error. */
/* ln(P) = -(Q/10)*ln(10) */
static const long double phred_multiplier = -0.2302585092994;

nuc_t decode_ref(char encoded);

/* Converts raw pileup read and qual strings into nuc_t's and log qualities */
int clean_fill(int read_depth, nuc_t ref, char* raw_reads, char* raw_quals, nuc_t* processed_reads, long double* processed_quals);

/*Based on maximum count, gets the alternate allele*/
nuc_t get_alt_allele(Locus* locus, nuc_t ref, int m);

/*Converts a nuc_t base to a human readable character */
char base2char(nuc_t base);


#endif
