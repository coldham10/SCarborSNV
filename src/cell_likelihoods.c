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
double l_P_amplification(nuc_t ref, nuc_t beta, int g, double l_P_err) {
    if (g == 0) {
        return (beta == ref) ? log(1 - exp(log(3) + log(l_P_err))) : l_P_err;
    }
    else if (g == 1) {
        if (beta == ref) {
            /* log((1-2*p_e)/2) */
            return log(1 - exp(log(2) + l_P_err)) - log(2);
        }
        else {
            /* (1/3)*(1- the above) */
            return (double)LSE2(-1 * log(6), l_P_err - log(3));
        }
    }
    else if (g == 2) {
        return (beta == ref) ? l_P_err : log(1-exp(l_P_err)) - log(3);
    }
    else {
        /*Error */
        return -1;
    }
}





/*Homozygous likelihood calculation */
long double homoz_cell_l(Cell_locus* cell, int g, nuc_t ref, double l_P_amp_err) {
    int k;
    nuc_t read;
    long double l_e;
    /* product found by summing in log space */
    long double product = 0;
    for (k = 0; k < cell->read_count; k++;) {
        read = cell->reads[k];
        l_e  = cell->quals[k];

        l_P_beta__g = l_P_amplification(ref, read, g);

    }

}

/*Heterozygous likelihood calculation */
long double hetz_cell_l(Cell_locus* cell, int g, nuc_t ref, double l_P_amp_err) {
    /*TODO*/
}

/*Given a Cell_locus, log-probability of amplification error and an 3-length array for results, compute P(g|d_ij)*/
int cell_likelihoods (Cell_locus* cell, long double* log_ls, nuc_t ref, double l_P_amp_err) {
    int i;
    if (cell->read_count < 1) {
        /*Indicate insufficient depth. XXX must check log_likelihoods to be positive*/
        log_ls[0] = log_ls[1] = log_ls[2] = -1;
        return 1;
    }
    /* g = 0,1,2*/
    for (i = 0; i < 3; i++) {
        if ((i % 2) == 0) {
            log_ls[i] = homoz_cell_l(cell, i, ref, l_P_amp_err);
        }
        else {
            log_ls[i] =  hetz_cell_l(cell, i, ref, l_P_amp_err);
        }
    }
    return 0;
}
