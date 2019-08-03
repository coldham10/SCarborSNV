#ifndef LOCUS_LIKELIHOODS_H_
#define LOCUS_LIKELIHOODS_H_

#include "sigma_priors.h"
#include "pileup_reader.h"

/*From reads stored in a Locus struct, fills in 2m+1 sigma likelihoods P(D|sig)
 * First cleans input of useless reads then uses DP to calculate sigma likelihoods
 * from cell genotype likelihoods */
int locus_likelihoods(Locus* locus, long double* locus_ls, prior_params_t* p);

#endif
