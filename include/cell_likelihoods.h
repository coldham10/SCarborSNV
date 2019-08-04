#ifndef CELL_LIKELIHOODS_H_
#define CELL_LIKELIHOODS_H_

/*Removes unknown and low quality cells then computes individual 
 * cell genotype log-likelihoods from read information and prior parameters.
 * A cell with no valid reads will get likelihoods of NAN. */
int cell_likelihoods (Cell_locus* cell, long double* log_ls, nuc_t ref, long double l_P_amp_err, long double l_P_ADO);

#endif