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
#ifndef POSTERIORS_H_
#define POSTERIORS_H_

#include "scarborsnv.h"

typedef struct {
    long double P_0, P_SNV;
    unsigned long pos;
    int m;
    char* seq_name;
    char* valid_cells;
    long double* simple_posteriors;
    long double*  phylo_posteriors;
    nuc_t ref, alt;
} Candidate;


/*Posterior distribution over sigma P(sigma|D) for a locus 
 * priors are original sigma priors and likelihoods  are locus likelihoods P(D|sigma)*/
int sigma_posteriors(long double* posteriors, long double* priors, long double* likelihoods, int m);
/*Given the posterior sigma distribution and cell likelihoods return posterior cell genotype probabilities*/
int cell_posteriors(long double* result, long double* sigma_dist, long double* cell_ls, int m);

#endif
