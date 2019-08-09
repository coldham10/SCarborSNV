#ifndef POSTERIORS_H_
#define POSTERIORS_H_

#include "scarborsnv.h"

typedef struct {
    long double P_0, P_SNV;
    unsigned long pos;
    int m;
    char* seq_name;
    char* valid_cells;
    long double* simple_posteriors;
    long double*  phylo_posteriors;
    nuc_t ref, alt;
} Candidate;


/*Posterior distribution over sigma P(sigma|D) for a locus 
 * priors are original sigma priors and likelihoods  are locus likelihoods P(D|sigma)*/
int sigma_posteriors(long double* posteriors, long double* priors, long double* likelihoods, int m);
/*Given the posterior sigma distribution and cell likelihoods return posterior cell genotype probabilities*/
int cell_posteriors(long double* result, long double* sigma_dist, long double* cell_ls, int m);

#endif
