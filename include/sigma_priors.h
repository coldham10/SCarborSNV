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
#ifndef SIGMA_PRIORS_H_
#define SIGMA_PRIORS_H_

typedef struct {
    long double l_P_tree; /*Calculated in sigma_priors.c */
    double c_thresh; /*Log of threshold for sigma=0 posterior to call candidates*/
    double l_lambda; /*Somatic mutation frequency*/
    double l_mu; /*Germline mutation frequency*/
    double l_P_H; /*Haploid event frequency. 0.09? this should be from empirical */
    /*FIXME get initial values for these two*/
    /*TODO add getopt options to change these following two*/
    double l_P_amp_err;
    double l_P_ADO;
    double l_P_clonal; /*0.51(?) refers to proportion of somatic genetic mutations that are public */
    int m;
} prior_params_t;

/* Returns the log prior probabilities for alternate allele count (sigma)
 * across a collection of m cells. Returned vector contains 2m+1 values
 * since sigma can range between 0 and 2m inclusive.
 */
int log_sigma_priors(prior_params_t* p, long double* l_priors);



#endif /*SIGMA_PRIORS_H_ */
