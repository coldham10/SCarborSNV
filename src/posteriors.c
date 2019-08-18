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
#include <stdlib.h>
#include "posteriors.h"
#include "math_utils.h"

int sigma_posteriors(long double* posteriors, long double* priors, long double* likelihoods, int m) {
    int i;
    long double sum;
    long double* products = malloc((2*m + 1) * sizeof(long double));
    for (i = 0; i < 2*m + 1; i++) {
        products[i] = priors[i] + likelihoods[i];
    }
    sum = LSE(products, 2*m + 1);
    for (i = 0; i < 2*m + 1; i++) {
        posteriors[i] = products[i] - sum;
    }
    free(products);
    return 0;
}

int cell_posteriors(long double* result, long double* sigma_dist, long double* cell_ls, int m) {
    /*TODO for improved efficiency, could only calculate this for non-blank cells?*/
    int s, cell_i;
    long double term_sum;
    long double terms[3];
    long double* to_sum0 = malloc((2*m + 1) * sizeof(long double));
    long double* to_sum1 = malloc((2*m + 1) * sizeof(long double));
    long double* to_sum2 = malloc((2*m + 1) * sizeof(long double));
    for (cell_i = 0; cell_i < m; cell_i++) {
        for (s = 0; s < 2*m + 1; s++) {
            /*Normalized P(d | sigma')*P(sigma')*/
            terms[0] = cell_ls[3*cell_i + 0] + l_binom_dist2(m, 0, s);
            terms[1] = cell_ls[3*cell_i + 1] + l_binom_dist2(m, 1, s);
            terms[2] = cell_ls[3*cell_i + 2] + l_binom_dist2(m, 2, s);
            term_sum = LSE(terms, 3);
            to_sum0[s] = terms[0] - term_sum + sigma_dist[s];
            to_sum1[s] = terms[1] - term_sum + sigma_dist[s];
            to_sum2[s] = terms[2] - term_sum + sigma_dist[s];
        }
        result[3*cell_i + 0] = LSE(to_sum0, 2*m + 1);
        result[3*cell_i + 1] = LSE(to_sum1, 2*m + 1);
        result[3*cell_i + 2] = LSE(to_sum2, 2*m + 1);
    }
    free(to_sum0); free(to_sum1); free(to_sum2); 
    return 0;
}
