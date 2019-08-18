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
#ifndef PILEUP_READER_H_
#define PILEUP_READER_H_

#include <stdio.h>
#include "scarborsnv.h"

typedef struct {
    nuc_t* reads;
    long double* quals;
    int read_count;
    int cell_position;
} Cell_locus;

typedef struct {
    char* sequence;
    Cell_locus* cells;
    unsigned long position;
    nuc_t ref_base;
} Locus;

/*Read n loci from instream into the array starting at loci.
 * returns the number of loci successfully allocated and read*/
int read_batch_loci(FILE* instream, Locus* loci, int n, int m);

/*Deletes anything alloc'd by read_batch_loci. WARNING: still need to 
 * free the array of Locus separately*/
int delete_locus_contents(Locus* loci, int n, int m);

#endif /*define PILEUP_READER_H_ */
