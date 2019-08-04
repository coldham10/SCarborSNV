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
