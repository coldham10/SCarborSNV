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
#ifndef LIKELIHOODS_H_
#define LIKELIHOODS_H_

#include "sigma_priors.h"

/*Removes unknown and low quality cells then computes individual 
 * cell genotype log-likelihoods from read information and prior parameters.
 * A cell with no valid reads will get likelihoods of NAN. */
int cell_likelihoods (Cell_locus* cell, long double* log_ls, nuc_t ref, long double l_P_amp_err, long double l_P_ADO);


/*Reads from cell likelihoods (in triplets) to calculate p(D|sig)
 * input is a 3*m array, output is size 2m+1 */
int locus_likelihoods(long double* cell_ls, long double* locus_ls, prior_params_t* p);

#endif
