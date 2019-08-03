#include <math.h>
#include "scarborsnv.h"
#include "math_utils.h"
#include "pileup_reader.h"

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
    /*FIXME needs l_P_ADO & l_P_amp_err. What initial values? should be in p0? check if monovar gives a p_e for table 6*/
    long double term1, term2;
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
