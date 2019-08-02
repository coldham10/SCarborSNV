#include "scarborsnv.h"
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

/*Homozygous likelihood calculation */
long double homoz_cell_l(Cell_locus* cell, nuc_t ref) {
    int k;
    /* product found by summing in log space */
    long double product = 0;
    for (k = 0; k < cell->read_count; k++;) {
        //
    }

}

/*Heterozygous likelihood calculation */
long double hetz_cell_l(Cell_locus* cell, nuc_t ref) {
    /*TODO*/
}

/*Given a Cell_locus and an 3-length array for results, compute P(g|d_ij)*/
int cell_likelihoods (Cell_locus* cell, long double* log_ls, nuc_t ref) {
    int i;
    if (cell->read_count < 1) {
        /*Indicate insufficient depth. XXX must check log_likelihoods to be positive*/
        log_ls[0] = log_ls[1] = log_ls[2] = -1;
        return 1;
    }
    /* g = 0,1,2*/
    for (i = 0; i < 3; i++) {
        if ((i % 2) == 0) {
            log_ls[i] = homoz_cell_l(cell, i, ref);
        }
        else {
            log_ls[i] = hetz_cell_l(cell, i, ref);
        }
    }
    return 0;
}
