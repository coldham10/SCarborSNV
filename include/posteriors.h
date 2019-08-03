#ifndef POSTERIORS_H_
#define POSTERIORS_H_

/*Posterior distribution over sigma P(sigma|D) for a locus */
int sigma_posteriors(long double* posteriors, long double* priors, long double* likelihoods, int m);

#endif
