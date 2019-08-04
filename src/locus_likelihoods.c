#include <stdlib.h>
#include <math.h>
#include "locus_likelihoods.h"
#include "cell_likelihoods.h"
#include "math_utils.h"


/*************************************************************
 * For a locus, the likelihood of sigma marginal on all cells,
 * P(D_i|sig), can be found using a dynamic programming 
 * technique similar to that of Monovar and Sciphi 
 * *********************************************************/

int locus_DP(long double* locus_ls, long double* cell_ls, prior_params_t* p) {
    int k, l;
    int m = p->m;
    /*Three terms of the equation to sum together*/
    long double* terms = malloc(3 * sizeof(long double));
    long double** DP_matrix = malloc(m * sizeof(long double*));
    /*Note k value here is one less than in paper*/
    for (k = 0; k < m; k++) {
        DP_matrix[k] = malloc((2*(k+1) + 1) * sizeof(long double));
    }
    /*Initial values*/
    DP_matrix[0][0] = cell_ls[0];
    DP_matrix[0][1] = logl(2) + cell_ls[1];
    DP_matrix[0][2] = cell_ls[2];
    /*Start DP algorithm*/
    for (k = 1; k < m; k++) {
        for (l = 0; l < 2*(k+1) + 1; l++) {
            if (l < 2*k + 1) {
                terms[0] = DP_matrix[k-1][l] + cell_ls[3*k];
            }
            else {
                terms[0] = logl(0);
            }
            if (l >= 1 && l < 2*k + 2) {
                terms[1] = logl(2) + DP_matrix[k-1][l-1] + cell_ls[3*k + 1];
            }
            else {
                terms[1] = logl(0);
            }
            if (l >= 2 && l < 2*k + 3) {
                terms[2] = DP_matrix[k][l-2] + cell_ls[3*k + 2];
            }
            else {
                terms[2] = logl(0);
            }

           DP_matrix[k][l] = LSE(terms, 3);
        }
    }
    /*Rescale for likelihood values*/
    for (l = 0; l <= 2*m; l++) {
        locus_ls[l] = DP_matrix[m-1][l] - log_binom(2*m, l);
    }
    /*Free matrix*/
    for (k = 0; k < m; k++) {
        free(DP_matrix[k]);
    }
    free(DP_matrix);
    free(terms);
    return 0;
}


int locus_likelihoods(long double* cell_ls, long double* locus_ls, prior_params_t* p) {
    int i;
    int n_valid_cells = p->m;
    for (i = 0; i < p->m; i++) {
        if (isnan(cell_ls[3*i])) {
            /*Neither allele was amplified, so likelihood set to p_ADO^2*/
            /*XXX is this legitamite */
            /*TODO test this with toy examples to see posterior sigma distribution vs prior */
            cell_ls[3*i] = cell_ls[3*i +1] = cell_ls[3*i +2] = 2 * p->l_P_ADO;
            n_valid_cells -= 1;
            continue;
        }
    }
    if (n_valid_cells == 0) {
        /*All likelihoods are just p_ADO^2m*/
        for (i = 0; i < p->m; i++) {
            locus_ls[i] = p->m * cell_ls[0];
        }
    }
    else {
        locus_DP(locus_ls, cell_ls, p);
    }
    return 0;
}
