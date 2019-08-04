#ifndef POSTERIORS_H_
#define POSTERIORS_H_

/*Posterior distribution over sigma P(sigma|D) for a locus 
 * priors are original sigma priors and likelihoods  are locus likelihoods P(D|sigma)*/
int sigma_posteriors(long double* posteriors, long double* priors, long double* likelihoods, int m);
/*Given the posterior sigma distribution and cell likelihoods return posterior cell genotype probabilities*/
int cell_posteriors(long double* result, long double* sigma_dist, long double* cell_ls, int m);

#endif
