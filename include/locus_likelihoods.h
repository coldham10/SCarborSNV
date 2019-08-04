#ifndef LOCUS_LIKELIHOODS_H_
#define LOCUS_LIKELIHOODS_H_

#include "sigma_priors.h"
#include "pileup_reader.h"

/*Reads from cell likelihoods (in triplets) to calculate p(D|sig)
 * input is a 3*m array, output is size 2m+1 */
int locus_likelihoods(long double* cell_ls, long double* locus_ls, prior_params_t* p);

#endif
