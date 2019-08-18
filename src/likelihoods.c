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
#include <math.h>
#include "scarborsnv.h"
#include "math_utils.h"
#include "pileup_reader.h"
#include "likelihoods.h"

/*************************************************************
 *
 * CELL LIKELIHOODS
 *
 ************************************************************/

/*Remove 'N' nucleotides, deletions, too low quality of reads*/
int prepare_reads(Cell_locus* cell) {
    int i;
    int n_removed = 0;
    if (cell->read_count < 1) {
        return 0;
    }
    for (i = 0; i < cell->read_count; i++) {
        /* This includes very low read quality, invalid, deleted, unknown, etc... */
        if ((cell->reads[i] >= INVALID_NUC) || (cell->quals[i] > logl(0.75)) ){
            n_removed += 1;
            continue;
        }
        cell->reads[i-n_removed] = cell->reads[i];
        cell->quals[i-n_removed] = cell->quals[i];
    }
    if (n_removed >= 1) {
        /* Mark end of 'clean' reads */
        cell->reads[cell->read_count-n_removed] = INVALID_NUC;
        cell->read_count -= n_removed;
    }
    return 0;
}

/*Returns log probability of amplifying to beta given the genotype
 * Using distribution from Zafar et al (2016), supplementary table 6*/
long double l_P_amplification(nuc_t ref, nuc_t beta, int g, long double l_P_err) {
    /* Since l_P_err does not change these are memoized for efficiency */
    static long double l_P1 = NAN;
    static long double l_P2 = NAN;
    static long double l_P3 = NAN;
    static long double l_P4 = NAN;
    /*Initialize on first call */
    if (isnan(l_P1)) {
        l_P1 = logl(1 - expl(logl(3) + l_P_err));
    }
    if (isnan(l_P2)) {
        /* log((1-2*p_e)/2) */
        l_P2 = logl(1 - expl(logl(2) + l_P_err)) - logl(2);
    }
    if (isnan(l_P3)) {
        /* (1/3)*(1- l_P2) */
        l_P3 = LSE2(-1 * logl(6), l_P_err - logl(3));
    }
    if (isnan(l_P4)) {
        l_P4 = logl(1-expl(l_P_err)) - logl(3);
    }
    /* return appropriate value */
    if (g == 0) {
        return (beta == ref) ?  l_P1 : l_P_err;
    }
    else if (g == 1) {
        return (beta == ref) ? l_P2 : l_P3;
    }
    else if (g == 2) {
        return (beta == ref) ? l_P_err : l_P4;
    }
    else {
        /*Error */
        return NAN;
    }
}


/*Homozygous likelihood calculation */
long double simple_cell_l(Cell_locus* cell, int g, nuc_t ref, long double l_P_amp_err) {
    int k;
    nuc_t read;
    long double l_e, l_P_beta__g, term1, term2;
    /* product found by summing in log space */
    long double product = 0;
    for (k = 0; k < cell->read_count; k++) {
        read = cell->reads[k];
        l_e  = cell->quals[k];

        l_P_beta__g = l_P_amplification(ref, read, g, l_P_amp_err);
        term1 = logl(1-expl(l_e)) + l_P_beta__g;
        term2 =  l_e - logl(3) + logl(1-expl(l_P_beta__g));
        product += LSE2(term1, term2);
    }

    return product;
}

/*Given a Cell_locus, log-probability of amplification error and an 3-length array for results, compute P(d_ij|g)*/
int cell_likelihoods (Cell_locus* cell, long double* log_ls, nuc_t ref, long double l_P_amp_err, long double l_P_ADO) {
    long double term1, term2;
    prepare_reads(cell);
    if (cell->read_count < 1) {
        /*Indicate insufficient depth. */
        log_ls[0] = log_ls[1] = log_ls[2] = NAN;
        return 1;
    }
    /*Homozygous are g=0,2*/
    log_ls[0] = simple_cell_l(cell, 0, ref, l_P_amp_err);
    log_ls[2] = simple_cell_l(cell, 2, ref, l_P_amp_err);
    /*Heterozygous case must consider with and without ADO*/
    term1 = l_P_ADO + LSE2(log_ls[0], log_ls[2]) - logl(2);
    term2 = logl(1-expl(l_P_ADO)) + simple_cell_l(cell, 1, ref, l_P_amp_err);
    log_ls[1] = LSE2(term1, term2);

    return 0;
}

/*************************************************************
 *
 * LOCUS LIKELIHOODS
 *
 ************************************************************/

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
        locus_ls[l] = DP_matrix[m-1][l] - l_binom_coeff(2*m, l);
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
