#ifndef CELL_LIKELIHOODS_H_
#define CELL_LIKELIHOODS_H_

/*Removes 'N' nucleotides, deletions, too low quality of reads*/
int prepare_reads(Cell_locus* cell);

/*Computes individual cell genotype log-likelihoods from read information and prior parameters*/
int cell_likelihoods (Cell_locus* cell, long double* log_ls, nuc_t ref, long double l_P_amp_err, long double l_P_ADO);

#endif
